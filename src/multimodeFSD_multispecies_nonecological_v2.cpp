/**************************************************************************
 * multimodeFSD_nonecological.cpp
 *
 * Source code for the non-ecological model described in Ashby and 
 * Boots "Multi-mode fluctuating selection in host-parasite coevolution". 
 * The text file "multimodeFSD_nonecological.txt" records:
 *
 * Column 1: time vector
 * Columns 2-(N+1): Host genotype abundances at each time point
 * Columns (N+2)-(2*N+1): Parasite genotype abundances at each time point
 *
 * The genotype ordering is output during the simulation
 *
 * The authors accept no responsibility or liability for any loss or damage
 * arising from the use of this code.
 *
 * Original code: Ben Ashby
 * 18/12/2016
 * v1.0 
 * Extension to a two-species model: Hanna Maerkle
 *************************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <getopt.h>
#include <stdio.h>

/***********************************
 * Model parameters
 ***********************************/
//#define LOCI 3 /* Number of loci */
//#define N (1 << LOCI) /* Number of genotypes (2^LOCI) */
#define MAXTIME 1e5 /* Length of simulation */

/***********************************
 * Runge-Kutta parameters
 ***********************************/
#define MAXSTEPS 1e6 /* Maximum number of steps for ODE solver */
#define INTERVAL 1e5 /* Interval for checking if the system is close to equilibrium */
#define EQTOL 0.0 /* Tolerance for checking if the system is close to equilibrium */
#define EPS 1e-6 /* ODE solver tolerance */
#define TINY 1e-6 /* Constant value for solver */


/////////////////////////////////////////////////////////////
/* Set standard values if not provided as input parameters */
////////////////////////////////////////////////////////////

unsigned int seed=1;
double CH1 = 0.15; /* Host maximum cost of resistance parameter 1 */
double CP1 = 0.15; /* Pathogen maximum cost of infectivity  1 */
double CH2 = 3.0; /* Host cost parameter 2 (>0:accelerating, =0:linear, <0:decelerating) */
double CP2 = 3.0; /* Parasite cost parameter 2 (>0:accelerating, =0:linear, <0:decelerating) */
double MINFREQ = 1e-6;

double BETAH = 1.0; /* Impact of infection on host fitness */
double BETAP = 1.0; /* Impact of infection on parasite fitness */

double SIGMA = 0.85; /* Resistance parameter */
double PHI = 0.2; /* Proportion of host species H1 as seen by the pathogen */

int tempi = 0;
double tempf = 0;
std::string temps = "";

int blockP = -1;
// Position in the array to initialize to be susceptible in host 1 
int block1 = -1;
// Position in the array to initialize to be susceptible in host 2 
int block2 = -1;
std::string prefix = "Non_ecological";
std::string initcond = "d";
int NLOCI = 3;
int excheck = 0;


/*************************************
 * Function prototypes
 *************************************/
void my_rungkut(double *INIT_POPH1, double *INIT_POPH2, double *INIT_POPP, double *CH_H1, double *CH_H2,
 double *CP_P2, double **Q1, double **Q2, int NH1TYPES, int NH2TYPES, int NPTYPES, std::string dfile);
void rkqs(double *XH1, double *XH2, double *P, double *DXH1DT, double *DXH2DT, double *DPDT, double *h, double *hnext, double *XH1SCALE, 
double *XH2SCALE, double *PSCALE, double **Q1, double **Q2, double *CH_H1, double *CH_H2, double *CP_P2, int NH1TYPES, int NH2TYPES, int NPTYPES);
/// @brief Function to run 
/// @param X Array of all current host genotype frequencies.
/// @param Y Array of all current pathogen genotype frequencies.
/// @param DXDT The result of evaluating the ODEs at the current time for all host.
/// @param DYDT The result of evaluating the ODEs at the current time for all pathogen
/// @param xout The calculated host genotype frequencies based on the Cash-Karp parameters
/// @param yout The calculated pathogen genotype frequencies based on the Cash-Karp parameters
/// @param xerr An array storing the local truncation error estimate for each host ODE.
/// @param yerr An array storing the local truncation error estimate for each pathogen ODE. 
/// @param h The tried step size.
/// @param Q The interaction matrix.
/// @param CH The base line fitness for each host genotype.
/// @param CP The base line fitness for each pathogen genotype
/// @param NPTYPES The total number of host and pathogen genotypes
void rkck(double *XH1, double *XH2, double *P, double *DXH1DT, double *DXH2DT, double *DPDT, double *x1out, double *x2out, double *pout, double *x1err, double *x2err, double *perr, double h, double **Q1, double **Q2, double *CH_H1, double *CH_H2, double *CP_P2, int NH1TYPES, int NH2TYPES, int NPTYPES);
/// @brief Function to calculate the ODEs for given values
/// @param X Array of the current host genotype frequencies
/// @param Y Array of the current pathogen genotype frequencies
/// @param Q Interaction matrix between host and pathogen frequencies
/// @param CH Array of the baseline fitness for each host genotype.
/// @param CP Array of the baseline fitness for each pathogen genotype.
/// @param DXDT Results of evaluating the ODEs for each host genotypes.
/// @param DYDT Results of evaluating the ODEs for each pathogen genotype.
/// @param NPTYPES Total number of genotypes in the host an the pathogen.
void dynamic(double *XH1, double *XH2, double *P, double **Q1, double **Q2, double *CH_H1, double *CH_H2, double *CP_P2, double *DXH1DT, double *DXH2DT, double *DPDT, int NH1TYPES, int NH2TYPES, int NPTYPES);
//void dynamic(double *X, double *XH2, double *Y, double **Q, double **Q2, double *CH, double *CH_H2, double *CP, double *DXDT, double *DXH2DT, double *DYDT, int NH2TYPES, int NTYPES);
/// @brief Function to find the maximum of two double values.
/// @param  First double value to compare.
/// @param  Second double value to compare.
/// @return The value of the larger of the two doubles compared.
double FMAX(double, double);
/// @brief Function to find the minimum of two double values
/// @param  first double to compare
/// @param  second double to compare
/// @return Return the value of the smaller of the two doubles being compared
double FMIN(double, double);
/// @brief Function to initialize the host genotypes, pathogen genotypes, calculate their range, their basefitness and the interaction matrix
/// @param CH Array for the baseline fitness of host 1/
/// @param CP Array for the baseline fitness for each pathogen genotype.
/// @param Q Array of pointers for storing the interaction matrix.
/// @param RANGE Array for storing the range of each genotype.
/// @param NTYPES Total numer of genotypes in the simulation.
/// @param NLOCI Total number of loci in the simulation.
/// @param geHfile File to which to write the information on all host genotypes.
/// @param gePfile File to which to write to the information on all pathogen genotypes.
void array_values(double *CH_H1, double *CH_H2, double *CP_P2, double **Q1, double **Q2, int *RANGEH1, int *RANGEH2, int *RANGEP, 
                    int NH1TYPES, int NH2TYPES, int NPTYPES, int NLOCI, int fscaleH1, int fscaleH2, int fscaleP, std::string geHfile, std::string geH2file, std::string gePfile);
/// @brief Help for the program
void PrintHelp();
/// @brief Function to read the command line arguments
/// @param argc 
/// @param argv 
void ProcessArgs(int argc, char** argv);
/// @brief Function to write all function parameters to an output file
/// @param CH1 Maximum loss of fitness due to genotype for the host.
/// @param CH2 Parameter determining the shape of the fitness function in the host.
/// @param CP1 Maximum loss of fitness due to genotype for the pathogen.
/// @param CP2 Parameter determining the shape of the fitness function in the pathogen.
/// @param SIGMA Sigma.
/// @param BETAH beta for the host.
/// @param BETAP beta for the pathogen.
/// @param block1 which locus to make monomorphic in the host.
/// @param block2 which locus to make monomorphic in host 2.
/// @param prefix The prefix to be used for all output files.
/// @param seed The seed to use for the simulation.
/// @param PHI The value of phi to use for the simulation.
/// @param paramfile The name of the output file for the parameters.
void writeparamsout (double CH1, double CH2, double CP1, double CP2, double SIGMA, double BETAH, double BETAP, int block1, int block2, std::string prefix, int seed, double PHI, std::string paramfile);
/// @brief Function to generate the dynamics to an output file (note happens continously during the simulation)
/// @param t Current time.
/// @param X Array of all current host genotype frequencies.
/// @param Y Array of all current pathogen gentype frequencies.
/// @param NPTYPES Total number of genotypes in the simulation.
/// @param out Name of the output file.
void write_dynamics(double t, double *XH1, double *XH2, double *P, int NH1TYPES, int NH2TYPES, int NPTYPES, std::ofstream& out);
/// @brief Function to make sure that none of the frequencies has gone negative due to numerical issues.
/// @param X Array of all current host frequencies.
/// @param Y Array of all current pathogen frequencies.
/// @param NPTYPES Total number of genotypes in the simulation.
void negativcheck(double *FREQS, int NTYPES);
/// @brief Function to make sure that all frequencies sum up to 1.
/// @param X Array of all current host genotype frequencies.
/// @param Y Array of all current pathogen gentoype frequnecies.
/// @param NTYPES Total number of genotypes in the simulation
void correctrounding(double *FREQS, int NTYPES, std::ofstream& out);
/// @brief Function to inialize the populations. 
/// @param INIT_FREQ Array in which the calculated initial frequencies will be stored.
/// @param TYPE_RANGE The range of each genotype.
/// @param NPTYPES The number of genotypes.
/// @param opt Indicator which initialization function should be used. Options are: f=Fixed initial frequency based on the range, r=random initial frequency, d=default as in Ashby and Boots.
void initialize_pop (double *INIT_FREQ, int *TYPE_RANGE, int NTYPES, int start,  char  opt);
/// @brief Function to print all genotypes to standard out
/// @param GENOTYPES Pointer array for the genotypes.
/// @param NPTYPES Number of genotypes in the simulation.
/// @param NLOCI Number of loci in the simulation.
void printgenotypes (int **GENOTYPES, int NTYPES, int NLOCI);
/// @brief Function to inialize all pathogen genotypes
/// @param GENOS the output.
/// @param NTYPES The number of genotypes to create.
/// @param NLOCI The number of loci.
void initialize_gts (int **GENOS, int NTYPES, int NLOCI, int block);

/// @brief Function to calculate the number of effective loci
/// @param NTYPES Total number of genotypes
/// @param NLOCI Total number of loci
/// @param GENOTYPES Genotype array
/// @param RANGE Return array
void calculateeffective (int NTYPES, int NLOCI, int **GENOTYPES, int *RANGE);
/// @brief Function to generate single element of the interaction matrix
/// @param GENOTYPE_HI Genotype of the host
/// @param GENOTYPE_PJ Genotype of the pathogen
/// @param SIGMA Sigma coefficient
/// @param NLOCI Number of loci in the simulation.
/// @return 
double calculate_Q(int *GENOTYPE_HI, int *GENOTYPE_PJ, double SIGMA, int NLOCI);
/// @brief Funtion to calculate the baseline fitness for a genotype in the absence of infection.
/// @param c_1 The maximum loss in fitness due to a certain genotype.
/// @param c_2 The shape parameter of the fitnes function.
/// @param normrange The range of the genotype.
/// @return 
double basefitness (double c_1, double c_2, double normrange);
/// @brief Function to write all host genotypes to a file at the beginning of the simulation.
/// @param GENOSH1 All genotypes in host species 1.
/// @param GENOSH2 All genotypes in host species 2.
/// @param GENOS All genotypes in the pathogen species.
/// @param RANGEH1 Range of each genotype in species 1.
/// @param RANGEH2 Range of each genotype in host species 2.
/// @param CH_H1 Baseline fitness of each host genotype in species 1.
/// @param CH_H2 Baseline fitness of each host genotypes in species 2.
/// @param Q1 Interaction matrix between all genotypes in host species 1 and all pathogen genotypes.
/// @param Q2 Interaction matrix between all genotypes in host species 2 and all pathogen genotypes.
/// @param NTYPES The number of genotypes in the simulation.
/// @param NLOCI The number of loci in the simulation.
/// @param geHfile The name of the output file.
void writehostgenotypefile (int **GENOSH1, int **GENOSH2, int **GENOSP, int *RANGEH1, int *RANGEH2, double *CH_H1, double *CH_H2, double **Q1, double **Q2, int NH1TYPES, int NH2TYPES, int NPTYPES, int NLOCI, int fscaleH1, int fscaleH2, std::string geHfile);
/// @brief Function to output the pathogen genotypes
/// @param GENOSH1 The genotypes of the first host species.
/// @param GENOSH2 The genotypes of the second host species.
/// @param GENOS The pathogen genotypes.
/// @param RANGEP The range of each pathogen genotype.
/// @param CP The basefitness of each pathogen genotype.
/// @param Q1 The interaction matrix with host species one.
/// @param Q2 The interaction matrix with host species two.
/// @param gePfile The name of the output file for the pathogen genotypes.
void writepathogengenotypefile(int **GENOSH1, int **GENOSH2, int **GENOSP, int *RANGEP, double *CP_P2, double **Q1, double **Q2, int NH1TYPES, int NH2TYPES, int NPTYPES, int NLOCI, std::string gePfile);
/// @brief Function to check if we are close to the equilibrium frequency (currently not used)
/// @param XMIN An array for the minimum frequency for a given host genotype across the two checkpoints. Note at the end of the check this will be set equal to the current frequency of the genotype.
/// @param YMIN An array for the minimum frequency for a given pathogen genotype across the two checkpoints. Note at the end of the check this will be set equal to the current frequency of the genotype.
/// @param XMAX An array for the maximum frequency for a given host genotype across the two checkpoints. Note at the end of the check this will be set equal to the current frequency of the genotype.
/// @param YMAX An array for the maximum frequency for a given pathogen genotype across the two checkpoints. Note at the end of the check this will be set equal to the current frequency of the genotype.
/// @param X The array of the current host genotype frequencies.
/// @param Y The array of the current pathogen genotype frequencoes.
/// @param t The current time
/// @param nextcheck The time when the next check for the equilibrium frequency should happen. If this checkpoint is not reached yet, nothing else is done
/// @param exitflag Are we done with the simulation.
/// @param NTYPES The number of genotypes in the simulation.
void eqcheck(double *XMIN, double *YMIN, double *XMAX, double *YMAX, double *X, double *Y, double *t, double *nextcheck, int *exitflag, int NTYPES);
void extinct_check(double *FREQS, int NTYPES, double minfreq);

int main (int argc, char* argv[]) {
    /* output files */
    std::string dfile = "";
    std::string geHfile = "";
    std::string geH2file = "";
    std::string gePfile = "";
    std::string paramfile = "";

    ProcessArgs(argc, argv);

    srand(seed);

    bool bpl = (blockP >= 0 );
    bool b1l = (block1 >= 0 );
    bool b2l = (block2 >= 0 );

    int NPTYPES = (bpl ? (1 << (NLOCI - 1)) : (1 << NLOCI)); 
    //std::cout << "NPTYPES " << NPTYPES << std::endl;
    //int NH1TYPES = (1 << NLOCI);
    int NH1TYPES = ( b1l ? (1 << (NLOCI - 1)) : (1 << NLOCI)); 
    //std::cout << "NH1TYPES " << NH1TYPES << std::endl;
    int NH2TYPES = ( b2l ? (1 << (NLOCI - 1)) : (1 << NLOCI));  // /* Number of host genotypes */
    //std::cout << "NH2TYPES " << NH2TYPES << std::endl;
    //std::cout << "bl1 " << int(b1l) << std::endl;
    //double INIT_POP[2*NTYPES], CH[NTYPES], CP[NTYPES];
    double INIT_POPH1[NH1TYPES], CH_H1[NH1TYPES];
    double INIT_POPH2[NH2TYPES], CH_H2[NH2TYPES];
    double INIT_POPP[NPTYPES], CP_P2[NPTYPES];
    int i;
    //int RANGE[NTYPES];
    int RANGEP[NPTYPES];
    int RANGEH1[NH1TYPES];
    int RANGEH2[NH2TYPES];


    // Allocate memory for an array of pointers
    //double** Q = new double*[NTYPES];
    double** Q1 = new double*[NH1TYPES];
    double** Q2 = new double*[NH2TYPES];

    // Allocate memory for each row
    /*
    for (int i = 0; i < NTYPES; ++i) {
        Q[i] = new double[NTYPES];
    }
    */

    for (i = 0; i < NH1TYPES; ++i) {
        Q1[i] = new double[NPTYPES];
    }

     // Allocate memory for each row
    for (i = 0; i < NH2TYPES; ++i) {
        Q2[i] = new double[NPTYPES];
    }   


    /* Print the read parameters to screen during run time */
    printf ("cH_1 = %f, cH_2 = %f, cP_1 = %f, cP_2= %f, sigma=%f, betaH=%f, betaP=%f, block1=%d, block2=%d, prefix=%s, seed=%d, phi= %f\n, initcond= %s\n, NLOCI= %d\n",
          CH1, CH2, CP1,CP2, SIGMA, BETAH, BETAP, block1, block2, prefix.c_str(),seed, PHI, initcond.c_str(), NLOCI);

    // Prepare the filenames 
    dfile = prefix + "_dynamics.txt"; // File to save the dynamics output to
    geHfile = prefix + "_genotypeH.txt"; // Host genotypes
    geH2file = prefix + "_genotypeH2.txt"; // Host genotypes
    gePfile = prefix + "_genotypeP.txt"; // Pathogen genotypes
    paramfile = prefix + "_params.txt"; // Parameters used in simulation

	// Store the input parameters in a file
	writeparamsout (CH1, CH2, CP1, CP2, SIGMA, BETAH, BETAP, block1, block2, prefix, seed, PHI, paramfile);
    
    /* Set constant array values */
    array_values(CH_H1, CH_H2, CP_P2, Q1, Q2, RANGEH1, RANGEH2, RANGEP, NH1TYPES, NH2TYPES, NPTYPES, NLOCI, NLOCI - int(b1l), NLOCI - int(b2l), NLOCI - int(bpl), geHfile, geH2file, gePfile);
    //exit(0);
    //initialize_pop(INIT_POP, RANGE, NTYPES, 0,  initcond[0]);
    //initialize_pop(INIT_POP, RANGE, NTYPES, NTYPES,  initcond[0]);
    initialize_pop(INIT_POPH1, RANGEH1, NH1TYPES, 0,  initcond[0]);
    initialize_pop(INIT_POPH2, RANGEH2, NH2TYPES, 0,  initcond[0]);
    initialize_pop(INIT_POPP, RANGEP, NPTYPES, 0,  initcond[0]);

    /*for(i=0; i < NH1TYPES; i++){
        std::cout << INIT_POPH1[i] << std::endl;
    }
    */

    /* Main routine */
    my_rungkut(INIT_POPH1, INIT_POPH2, INIT_POPP, CH_H1, CH_H2, CP_P2, Q1, Q2, NH1TYPES, NH2TYPES, NPTYPES, dfile);
    
    return 0;
}


/*****************************************
 * Print the help
 ****************************************/
void PrintHelp()
{
    std::cout <<
            "--cH_1 <double>:			Set cH_1\n"
            "--cH_2 <double>:			Set cH_2\n"
            "--cP_1 <double>:			Set cP_1\n"
            "--cP_2 <double>:			Set cP_2\n"
            "--fprefix <string>:		Set prefix for all output files\n"
            "--block 1 <int>:			Which position to block H1\n"
            "--block 2 <int>:			Which position to block H2\n"
            "--sigma <double>:			Specify sigma\n"
            "--betaH <double>:			set betaH\n"
            "--betaP <double>:			set betaP\n"
            "--phi <phi>:				set phi\n"
            "--init <string>:			set init\n"
            "--seed <int>:				set seed\n"
            "--nloci <int>:             number of loci in the simulation\n"
            "--help:					Show help\n";
    exit(1);
}


/*****************************************
 * Command line argument parser
 ****************************************/
void ProcessArgs(int argc, char** argv)
{
    const char* const short_opts = "g:i:p:q:f:b:c:s:d:e:a:r:z:n:y:h";
    const option long_opts[] = {
            {"cH_1", required_argument, nullptr, 'g'},
            {"cH_2", required_argument, nullptr, 'i'},
            {"cP_1", required_argument, nullptr, 'p'},
            {"cP_2", required_argument, nullptr, 'q'},
            {"fprefix", required_argument, nullptr, 'f'},
            {"block1", required_argument, nullptr, 'b'},
            {"block2", required_argument, nullptr, 'c'},
            {"sigma", required_argument, nullptr, 's'},
            {"betaH" , required_argument, nullptr, 'd'},
            {"betaP",  required_argument, nullptr, 'e'},
            {"phi", required_argument, nullptr, 'a'},
            {"seed", required_argument, nullptr, 'r'},
            {"init", required_argument, nullptr, 'z'},
            {"nloci", required_argument, nullptr, 'n'},
            {"excheck", required_argument, nullptr, 'y'},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, 0, nullptr, 0} // End of options
    };
    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {
			case 'b':
                tempi = block1;
                block1 = std::atoi(optarg);
                if( (block1 < -1) | (block1 > 2) ){
                    fprintf (stderr, "Option block1 can be only in {-1,0,1,2}. Resetting value to default value %i\n", tempi);
                    block1 = tempi;
                    }
                break;       
            case 'c':
                tempi = block2;
                block2 = std::atoi(optarg);
                if( (block2 < -1) | (block2 > 2) ){
                    fprintf (stderr, "Option block2 can be only in {-1,0,1,2}. Resetting value to default value %i\n", tempi);
                    block2 = tempi;
                    }
                break;            

            case 'g':
                tempf = CH1;
                CH1 = std::atof(optarg);
                if( (CH1 < 0) | (CH1 > 1) ){
                    fprintf (stderr, "Option cH_1 needs to be in [0,1]. Resetting value to default value %f\n", tempf);
                    CH1 = tempf;
                    }
                break;

            case 'i':
                tempf = CH2;
                CH2 = std::atof(optarg);
                if( (CH2 < (-3)) | (CH2 > 3) ){
                    fprintf (stderr, "Option cH_2 needs to be in [-3,3]. Resetting value to default value %f\n", tempf);
                    CH2 = tempf;
                    }
                break;
            case 'p':
                tempf = CP1;
                CP1 = std::atof(optarg);
                if( (CP1 < 0) | (CP1 > 1) ){
                    fprintf (stderr, "Option cP_1 needs to be in [0,1]. Resetting value to default value %f\n", tempf);
                    CP1 = tempf;
                    }
                break;

            case 'q':
                tempf = CP2;
                CP2 = std::atof(optarg);
                if( (CP2 < (-3)) | (CP2 > 3) ){
                    fprintf (stderr, "Option cP_2 cH_2 needs to be in [-3,3]. Resetting value to default value %f\n", tempf);
                    CP2 = tempf;
                    }
                break;       
            case 's':
                tempf = SIGMA;
                SIGMA = std::atof(optarg);
                if( (SIGMA < 0) | (SIGMA > 1) ){
                    fprintf (stderr, "Option sigma needs to be in [0,1]. Resetting value to default value %f\n", tempf);
                    SIGMA = tempf;
                    }
                break;
            case 'd':
                tempf = BETAH;
                BETAH = std::atof(optarg);
                if( (BETAH < 0) | (BETAH > 1) ){
                    fprintf (stderr, "Option betaH needs to be in [0,1]. Resetting value to default value %f\n", tempf);
                    BETAH = tempf;
                    }
                break;           
            case 'e':
                tempf = BETAP;
                BETAP = std::atof(optarg);
                if( (BETAP < 0) | (BETAP > 1) ){
                    fprintf (stderr, "Option betaP needs to be in [0,1]. Resetting value to default value %f\n", tempf);
                    BETAP = tempf;
                    }
                break;                                              
            case 'a':
                tempf = PHI;
                PHI = std::atof(optarg);
                if( (PHI < 0) | (PHI > 1) ){
                    fprintf (stderr, "Option phi needs to be in [0,1]. Resetting value to default value %f\n", tempf);
                    PHI = tempf;
                    }
                break;   
            case 'f':
                prefix = optarg;
                break;
            case 'z':
                temps = initcond;
                initcond = optarg;
                if( (initcond != "f") && (initcond != "r") && (initcond != "d") ){
                    fprintf (stderr, "Unknown argument %s for option --init. Resetting to default value %s\n", initcond.c_str(),temps.c_str());
                    }                
                break;
            case 'n':
                tempi = NLOCI;
                NLOCI = std::atoi(optarg);
                if( (NLOCI < 0) | (NLOCI > 5) ){
                    fprintf (stderr, "Option --loci needs to be in [1..5]. Resetting value to default value %i\n", tempi);
                    NLOCI = tempi;
                    }
                break;               
            case 'r':
                seed = std::atoi(optarg);
                break;
            case 'y':
                excheck = std::atoi(optarg);
                break;                

        case 'h': // -h or --help
        case '?': // Unrecognized option
        default:
            PrintHelp();
            break;
        }
    }
}

/*****************************************
 * Keep a log of the parameters
 ****************************************/
void writeparamsout (double CH1, double CH2, double CP1, double CP2, double SIGMA, double BETAH, double BETAP, int block1, int block2, std::string prefix, int seed, double PHI, std::string paramfile){
	std::ofstream outpar(paramfile);
	
	outpar << "cH_1 " << CH1 << 
        "\ncH_2 " << CH2 << 
        "\ncP_1 " << CP1  << 
        "\ncP_2 " << CP2 << 
        "\nsigma " << SIGMA <<  
        "\nbetaH " << BETAH <<
        "\nbetaP " << BETAP << 
        "\nblock1 " << block1 << 
        "\nblock2 " << block2 <<
        "\nprefix " << prefix << 
        "\nseed " << seed << 
        "\nphi " << PHI << 
        "\ninitcond " << initcond 
        << std::endl;
}

/*****************************************
 
******************************************/

void writehostgenotypefile (int **GENOSH1, int **GENOSH2, int **GENOSP, int *RANGEH1, int *RANGEH2, double *CH_H1, double *CH_H2, double **Q1, double **Q2, int NH1TYPES, int NH2TYPES, int NPTYPES, int NLOCI, int fscaleH1, int fscaleH2, std::string geHfile){
	std::ofstream genooutH(geHfile);
	int i, j, k;
		
	genooutH << "Host " << "ID " << "Genotype";
    
    for(j=0; j < NPTYPES; j++){
            genooutH << " P_";
            for(k = 0; k < NLOCI; k++){
                genooutH << std::to_string(GENOSP[j][k]);
            }
    }    
    genooutH << " Basic_fitness " << "Range" << std::endl;

    for(i = 0; i < NH1TYPES; i++){
        genooutH << "H1 " << i << " ";
        for(k=0; k < NLOCI; k++){
            genooutH << std::to_string(GENOSH1[i][k]);
        }
        for(j = 0; j < NPTYPES; j++){
            genooutH << " " << Q1[i][j];
        }
        genooutH << " " << CH_H1[i] << " " << std::to_string((double)RANGEH1[i]/(double)(fscaleH1)) << std::endl; 
    }
    

    for(i = 0; i < NH2TYPES; i++){
        genooutH << "H2 " << i << " ";
        for(k=0; k < NLOCI; k++){
            genooutH << std::to_string(GENOSH2[i][k]);
        }
        for(j=0; j < NPTYPES; j++){
            genooutH << " " << Q2[i][j];
        } 
        genooutH << " " << CH_H2[i] << " " << std::to_string((double)RANGEH2[i]/(double)(fscaleH2)) << std::endl;
    }
}


void writepathogengenotypefile (int **GENOSH1, int **GENOSH2, int **GENOSP, int *RANGEP, double *CP, double **Q1, double **Q2, int NH1TYPES, int NH2TYPES, int NPTYPES, int NLOCI, std::string gePfile) {
	int i, k, j;
	
	std::ofstream genooutP(gePfile);


    genooutP << "Para " << "ID " << "Genotype";
    for(i=0; i < NH1TYPES; i++){
            genooutP << " H1_";
            for(k=0; k < NLOCI; k++){
                genooutP << std::to_string(GENOSH1[i][k]);
            }
    }  

    for(i=0; i < NH2TYPES; i++){
            genooutP << " H2_";
            for(k=0; k < NLOCI; k++){
                genooutP << std::to_string(GENOSH2[i][k]);
            }
    }      
    genooutP << " Basic_fitness " << "Range" << std::endl;

    for(j=0; j < NPTYPES; j++){
        genooutP << "P " << j << " ";
        for(k=0; k < NLOCI; k++){
            genooutP << std::to_string(GENOSP[j][k]);
        }
        for(i=0; i < NH1TYPES; i++){
            genooutP << " " << Q1[i][j];
        } 
        for(i=0; i < NH2TYPES; i++){
            genooutP << " " << Q2[i][j];
        } 

        genooutP << " " << CP[j] << " " << std::to_string((double)RANGEP[j]/(double)(NLOCI)) << std::endl;
    }
}
		

/*****************************************
 * Function to write the dynamics to file
 ****************************************/
void write_dynamics(double t, double *XH1, double *XH2,  double *P, int NH1TYPES, int NH2TYPES, int NPTYPES, std::ofstream& out){
	int i;
	
	/* Update output */
    out << t << " ";

    for (i=0; i < NH1TYPES; i++) {
        out << XH1[i] << " ";
    }

    for (i=0; i < NH2TYPES; i++) {
        out << XH2[i] << " ";
    }
    //std::cout << "sdfs" << std::endl;

    for (i=0; i < NPTYPES; i++) {
        out << P[i] << " ";
    }
    

    out << "\n";
}

/*****************************************
 * Make sure that none of the frequencies goes negative
 ****************************************/
void negativcheck(double *FREQS, int NTYPES){
	int i;
	
	for(i=0; i < NTYPES; i++)
    {
            FREQS[i] = FMAX(FREQS[i],0);
    }
}

/*****************************************
 * Make sure that the frequencies sum up to 1
 ****************************************/
void correctrounding(double *FREQS, int NTYPES, std::ofstream& out){
	double total = 0;
    int i;

    /* sum up host frequencies*/
    for (i=0; i < NTYPES; i++){
        total += FREQS[i];
    }

    if(total == 0){
        std::cout << "Stoped with total criterium " << std::endl;
        out.close();
        exit(2);
    }

    /* Renormalize the host*/
    for (i=0; i < NTYPES; i++){
        FREQS[i] /= total;
    }

}


/*****************************************
 * Initialize populations
 ****************************************/
void initialize_pop(double *INIT_FREQ, int *RANGE, int NTYPES, int start, char opt){
	int i;

	/* Set initial population */
    double total = 0;

    switch(opt) {
        case 'f':
            std::cout << "Running case: f" << std::endl;
            for(i = 0; i < NTYPES; i++) {
                INIT_FREQ[i + start] = pow(10,-(double)RANGE[i]);
                total += INIT_FREQ[i + start];
            }
            break;
        case 'r':
            std::cout << "Running case: r" << std::endl;
            for(i = 0; i < NTYPES; i++) {
                INIT_FREQ[i+ start] = FMAX(0.1,((double)rand()/RAND_MAX));
                total += INIT_FREQ[i + start];
            }
            break;
        case 'd':
            std::cout << "Running case: d" << std::endl;
            for(i = 0; i < NTYPES; i++) {
                INIT_FREQ[i + start] = FMAX(0.1,((double)rand()/RAND_MAX))*pow(10,-(double)RANGE[i]);
                total += INIT_FREQ[i + start];
            }
            break;
        default:
            std::cout << "Running default case" << std::endl;
            for(i = 0; i < NTYPES; i++) {
                INIT_FREQ[i + start] = FMAX(0.1,((double)rand()/RAND_MAX))*pow(10,-(double)RANGE[i]);
                total += INIT_FREQ[i + start];
            }
    }
    
    /* Normalise */
    for(i=0; i < NTYPES; i++) {
        INIT_FREQ[i + start] /= total;
    }
}

/*****************************************
 * Print genotypes
 ****************************************/
void printgenotypes (int **GENOTYPES, int NTYPES, int NLOCI){
	int i,j;
	std::string add=" ";
	
	for(i=0; i < NTYPES; i++){
		for(j=0; j < NLOCI; j++){
			std::cout << GENOTYPES[i][j] << add;
		}
		std::cout << "\n";
	}
}
	

/*****************************************
 * ODE solvre
 ****************************************/
void my_rungkut (double *INIT_POPH1, double *INIT_POPH2, double *INIT_POPP, double *CH_H1, double *CH_H2, double *CP, double **Q1, double **Q2, int NH1TYPES, int NH2TYPES, int NPTYPES, std::string dfile){
    
    //double X[NTYPES], DXDT[NTYPES],  XSCALE[NTYPES], XMIN[NTYPES], XMAX[NTYPES];
    double XH1[NH1TYPES], DXH1DT[NH1TYPES], XH1SCALE[NH1TYPES]; //, XMINH1[NH1TYPES], XMAXH1[NH1TYPES];
    double XH2[NH2TYPES], DXH2DT[NH2TYPES], XH2SCALE[NH2TYPES]; //, XMINH2[NH2TYPES], XMAXH2[NH2TYPES];
    //double Y[NTYPES], DYDT[NTYPES], YSCALE[NTYPES], YMIN[NTYPES], YMAX[NTYPES];
    double P[NPTYPES], DPDT[NPTYPES], PSCALE[NPTYPES]; //, PMIN[NPTYPES], PMAX[NPTYPES];

    double hnext[1], h[1];
    double t; //, nextcheck;
    int i, j, exitflag, count, maxsteps;
    //char filename[100];
    
    /* Create output file */
    std::cout << "Writing dynamics to dfile: " << dfile << "\n";

    std::ifstream inhile(dfile);
    if(inhile.good()) {
        std::cout << "File already exists, deleting ...\n";
        remove(dfile.c_str());
    }
    std::ofstream out(dfile, std::ios::app);
    if (!out){
        std::cout << "Cannot create file\n";
        exit(1);
    }
    
    /* Other parameters */
    exitflag = 1;
    count = 0;
    h[0] = 1e-3;
    hnext[0] = 1e-3;
    t = 0;
    //nextcheck = (int)INTERVAL;
    maxsteps = (int)MAXSTEPS;
    
    /* Initialise populations - start with all present */
    /*
    for(i=0; i < NTYPES; i++){
        X[i] = INIT_POP[i];
        XMIN[i] = X[i];
        XMAX[i] = X[i];
        Y[i] = INIT_POP[i+NTYPES];
        YMIN[i] = Y[i];
        YMAX[i] = Y[i];
    }
    */

    for(i=0; i < NH1TYPES; i++){
        XH1[i] = INIT_POPH1[i];
        //XMINH1[i] = XH1[i];
        //XMAXH1[i] = XH1[i];
    }

    for(i=0; i < NH2TYPES; i++){
        XH2[i] = INIT_POPH2[i];
        //XMINH2[i] = XH2[i];
        //XMAXH2[i] = XH2[i];
    }

    for(j=0; j < NPTYPES; j++){
        P[j] = INIT_POPP[j];
        //PMIN[j] = P[j];
        //PMAX[j] = P[j];
    }

    
    write_dynamics(t, XH1, XH2, P, NH1TYPES, NH2TYPES, NPTYPES, out);
    
    /* Main loop: */
    do{
        /* This ensures the final step lands us on the final time point */
        if( (1.1 * hnext[0]) > (MAXTIME-t)){
            hnext[0] = MAXTIME-t;
            h[0] = MAXTIME-t;
            t = MAXTIME;
            exitflag = 0;
        }
        else{
            h[0] = hnext[0];
            t += h[0];
        }
        if(t >= MAXTIME) {
            t = MAXTIME;
            exitflag = 0;
        }
        
        /* This is where your equations are first solved */
        dynamic(XH1, XH2, P, Q1, Q2, CH_H1, CH_H2, CP, DXH1DT, DXH2DT, DPDT, NH1TYPES, NH2TYPES, NPTYPES);
                
        /* Adjust the step size to maintain accuracy */
        /*
        for (i=0; i < NTYPES; i++){
            X[i] = FMAX(X[i],0);
            Y[i] = FMAX(Y[i],0);
            XSCALE[i] = fabs(X[i]) + fabs(DXDT[i] * (*h)) + TINY;
            YSCALE[i] = fabs(Y[i]) + fabs(DYDT[i] * (*h)) + TINY;
        } 
        */       

        for(i=0; i < NH1TYPES; i++){
            XH1[i] = FMAX(XH1[i],0);
            XH1SCALE[i] = fabs(XH1[i]) + fabs(DXH1DT[i] * (*h)) + TINY;
        }


        for(i=0; i < NH2TYPES; i++){
            XH2[i] = FMAX(XH2[i],0);
            XH2SCALE[i] = fabs(XH2[i]) + fabs(DXH2DT[i] * (*h)) + TINY;
        }

        for(i=0; i < NPTYPES; i++){
            P[i] = FMAX(P[i],0);
            PSCALE[i] = fabs(P[i]) + fabs(DPDT[i] * (*h)) + TINY;
        }

        rkqs(XH1, XH2, P, DXH1DT, DXH2DT, DPDT, h, hnext, XH1SCALE, XH2SCALE, PSCALE, Q1, Q2, CH_H1, CH_H2, CP, NH1TYPES, NH2TYPES, NPTYPES);
        
        if(excheck){
            extinct_check(XH1, NH1TYPES, MINFREQ);
            extinct_check(XH2, NH2TYPES, MINFREQ);
            extinct_check(P, NPTYPES, MINFREQ);
        }

        //negativcheck(X, NTYPES);
        //std::cout << "Negative check" << std::endl;
        negativcheck(XH1, NH1TYPES);
        negativcheck(XH2, NH2TYPES);
        //negativcheck(Y, NTYPES);
        negativcheck(P, NPTYPES);


        //correctrounding(X, NTYPES);
        correctrounding(XH1, NH1TYPES, out);
        correctrounding(XH2, NH2TYPES, out);
        //correctrounding(Y, NTYPES);
        correctrounding(P, NPTYPES, out);

        count++;
        
        //std::cout << t << std::endl;
        write_dynamics(t, XH1, XH2, P, NH1TYPES, NH2TYPES, NPTYPES, out);

        /* Check if we're close to equilibrium */
        //if(EQTOL > 0){
        //    eqcheck(XMIN, YMIN, XMAX, YMAX, X, Y, &t, &nextcheck, &exitflag, NTYPES);
        //}

    }while( count<(maxsteps-1) && t <= MAXTIME && exitflag);
}

void eqcheck(double *XMIN, double *YMIN, double *XMAX, double *YMAX, double *X, double *Y, double *t, double *nextcheck, int *exitflag, int NTYPES){
    int i;

    std::cout << "exitflag: " << *exitflag << std::endl;
    std::cout << "t: " << *t << std::endl;
    std::cout << "nextcheck: " << *nextcheck << std::endl;

    for (i=0; i < NTYPES; i++){
        XMIN[i] = FMIN(XMIN[i], X[i]);
        XMAX[i] = FMAX(XMAX[i], X[i]);
        YMIN[i] = FMIN(YMIN[i], Y[i]);
        YMAX[i] = FMAX(YMAX[i], Y[i]);
    }
            
    if( t > nextcheck){
        std::cout << "In next check condition" << std::endl;
        *exitflag = 0;
        for (i=0; i < NTYPES; i++){
            if(fabs(XMAX[i]-XMIN[i])>EQTOL || (YMAX[i]-YMIN[i]) > EQTOL){
                *exitflag = 1;
                std::cout << "OUTBREAK: " << i <<  std::endl;
                break;
            }
        }
    
        
        if(*exitflag == 0){
            *t = MAXTIME;
            std::cout << "Exit zero condition" << std::endl;
            return;
        }
                
        *nextcheck += INTERVAL;
        for (i=0; i < NTYPES; i++){
            XMIN[i] = X[i];
            XMAX[i] = X[i];
            YMIN[i] = Y[i];
            YMAX[i] = Y[i];
        }
    }
    
    //std::cout << "end exitflag: " << *exitflag << std::endl;
    //std::cout << "end t: " << *t << std::endl;
    //std::cout << "end nextcheck: " << *nextcheck << std::endl;   

}


/***********************************
 * Generate the adaptive step-size
 ***********************************/
void rkqs(double *XH1, double *XH2, double *P, double *DXH1DT, double *DXH2DT, double *DPDT, double *h, double *hnext, double *XH1SCALE, double *XH2SCALE, double *PSCALE, double **Q1, double **Q2, double *CH_H1, double *CH_H2, double *CP, int NH1TYPES, int NH2TYPES, int NPTYPES)
{
    //double xtemp[NTYPES], xerr[NTYPES];
    double x1temp[NH1TYPES], x1err[NH1TYPES];
    double x2temp[NH2TYPES], x2err[NH2TYPES];
    //double ytemp[NTYPES], yerr[NTYPES];
    double ptemp[NPTYPES], perr[NPTYPES];
    double htemp, errmax;
    int i; //, j;
    double phi1 = PHI;
    double phi2 = 1.0 - PHI;
    
    for(;;)
    {
        rkck(XH1, XH2, P, DXH1DT, DXH2DT, DPDT, x1temp, x2temp, ptemp, x1err, x2err, perr,  *h, Q1, Q2, CH_H1, CH_H2, CP, NH1TYPES, NH2TYPES, NPTYPES);        
        errmax = 0.0;
        /*
        for(i=0; i < NTYPES ;i++){
            errmax = FMAX(errmax, fabs(xerr[i]/XSCALE[i]));
            errmax = FMAX(errmax, fabs(yerr[i]/YSCALE[i]));
        }
        */
        for(i=0; i < NPTYPES ;i++){
            errmax = FMAX(errmax, fabs(perr[i]/PSCALE[i]));
        }
        
        if(phi1 > 0.0){
            for(i=0; i < NH1TYPES ;i++){
                errmax = FMAX(errmax, fabs(x1err[i]/XH1SCALE[i]));
            }
        }

        if(phi2 > 0.0){
            for(i=0; i < NH2TYPES ;i++){
                errmax = FMAX(errmax, fabs(x2err[i]/XH2SCALE[i]));
            }
        }

        errmax /= EPS;

        if(errmax<=1.0) break; /* We accept the current step size to be sufficient error wise and move forward */
        //std::cout << "Lost in loop" << std::endl;

        htemp= 0.9*(*h) * pow(errmax, -0.25);
        *h = (*h>=0.0 ? FMAX(htemp, 0.1 * (*h)) : FMIN(htemp, 0.1 * (*h)));
        //std::cout << "new size " << *h << std::endl;
        //exit(0);
    }

    if(errmax > 1.89E-4) {
        *hnext = 0.9* (*h) * pow(errmax, -0.2);
    }
    else {
        *hnext= 5.0 * (*h);
    }
    
    /*
    for(i=0; i < NTYPES; i++){
        X[i]= xtemp[i];
        Y[i]= ytemp[i];
    }
    
    std::cout << "MADE IT HERE" << std::endl;
    */

    for(i=0; i < NH1TYPES; i++){
        XH1[i] = x1temp[i];
    }

    for(i=0; i < NH2TYPES; i++){
        XH2[i] = x2temp[i];
    }

    for(i=0; i < NPTYPES; i++){
        P[i] = ptemp[i];
    }
}

/***********************************
 * Standard RK solver
 ***********************************/
void rkck(double *XH1, double *XH2, double *P, double *DXH1DT, double *DXH2DT, double *DPDT, double *x1out, double *x2out, double *pout, double *x1err, double *x2err, double *perr, double h, double **Q1, double **Q2, double *CH_H1, double *CH_H2, double *CP, int NH1TYPES, int NH2TYPES, int NPTYPES)
{
    int i; //, j;
    //double xk1[NTYPES], xk2[NTYPES], xk3[NTYPES], xk4[NTYPES], xk5[NTYPES], xk6[NTYPES], xtemp[NTYPES];
    //double x1k1[NH1TYPES], 
    double x1k2[NH1TYPES], x1k3[NH1TYPES], x1k4[NH1TYPES], x1k5[NH1TYPES], x1k6[NH1TYPES], x1temp[NH1TYPES];
    //double x2k1[NH2TYPES], x2k2[NH2TYPES], x2k3[NH2TYPES], x2k4[NH2TYPES], x2k5[NH2TYPES], x2k6[NH2TYPES], x2temp[NH2TYPES];
    double x2k2[NH2TYPES], x2k3[NH2TYPES], x2k4[NH2TYPES], x2k5[NH2TYPES], x2k6[NH2TYPES], x2temp[NH2TYPES];
    //double yk1[NTYPES], yk2[NTYPES], yk3[NTYPES], yk4[NTYPES], yk5[NTYPES], yk6[NTYPES], ytemp[NTYPES];
    //double pk1[NTYPES], pk2[NTYPES], pk3[NTYPES], pk4[NTYPES], pk5[NTYPES], pk6[NTYPES], ptemp[NTYPES];
    double pk2[NPTYPES], pk3[NPTYPES], pk4[NPTYPES], pk5[NPTYPES], pk6[NPTYPES], ptemp[NPTYPES];
    /* Kept for completeness
    const float a2 = 0.2;
    const double a3 = 0.3;
    const double a4 = 0.6;
    const double a5 = 1.0;
    const double a6 = 0.875;
    */
    const double b21 = 0.2;
    const double b31 = 3.0/40.0;
    const double b32 = 9.0/40.0;
    const double b41 = 0.3;
    const double b42 = -0.9;
    const double b43 = 1.2;
    const double b51 = -11.0/54.0;
    const double b52 = 2.5;
    const double b53 = -70.0/27.0;
    const double b54 = 35.0/27.0;
    const double b61 = 1631.0/55296.0;
    const double b62 = 175.0/512.0;
    const double b63 = 575.0/13824.0;
    const double b64 = 44275.0/110592.0;
    const double b65 = 253.0/4096.0;
    const double c1 = 37.0/378.0;
    //const double c2 = 0;
    const double c3 = 250.0/621.0;
    const double c4 = 125.0/594.0;
    const double c5 = 0.0;
    const double c6 = 512.0/1771.0;
    const double dc1 = c1 - 2825.0/27648.0;
    //const double dc2 = 0.0;
    const double dc3 = c3 - 18575.0/48384.0;
    const double dc4 = c4 - 13525.0/55296.0;
    const double dc5 = c5 - 277.00/14336.0;
    const double dc6 = c6 - 0.25;
    
    /*
    for(i=0;i < NTYPES; i++){
        xtemp[i] = X[i] + b21*h*DXDT[i];
        ytemp[i] = Y[i] + b21*h*DYDT[i];
    }
    */

    for(i=0; i < NH1TYPES; i++){
        x1temp[i] = XH1[i] + b21*h*DXH1DT[i];
        /*
        std::cout << "i: " << i << std::endl;
        std::cout << "DXH1DT " << DXH1DT[i] << std::endl;
        std::cout << "XH1 " << XH1[i] << std::endl;
        std::cout << "x1temp " << x1temp[i] << std::endl;
        */
    }

    for(i=0; i < NH2TYPES; i++){
        x2temp[i] = XH2[i] + b21*h*DXH2DT[i];
        /*
        std::cout << "i: " << i << std::endl;
        std::cout << "DXH2DT " << DXH2DT[i] << std::endl;
        std::cout << "XH2 " << XH2[i] << std::endl;
        std::cout  << "x2temp " << x2temp[i] << std::endl;
        */
    }


    for(i=0; i < NPTYPES; i++){
        ptemp[i] = P[i] + b21*h*DPDT[i];
        /*
        std::cout << "i: " << i << std::endl;
        std::cout << "DPDT " << DPDT[i] << std::endl;
        std::cout << "P " << P[i] << std::endl;
        std::cout  << "ptemp " << ptemp[i] << std::endl;
        */
    }

    dynamic(x1temp, x2temp, ptemp, Q1, Q2, CH_H1, CH_H2, CP, x1k2, x2k2, pk2, NH1TYPES, NH2TYPES, NPTYPES);
    //exit(0);
    /*   
    for(i=0; i < NTYPES; i++){
        xtemp[i] = X[i] + h * (b31*DXDT[i] + b32*xk2[i]);
        ytemp[i] = Y[i] + h * (b31*DYDT[i] + b32*yk2[i]);
    }
    */ 

    for(i=0; i < NH1TYPES; i++){
        x1temp[i] = XH1[i] + h * (b31*DXH1DT[i] + b32*x1k2[i]);
    }

    for(i=0; i < NH2TYPES; i++){
        x2temp[i] = XH2[i] + h * (b31*DXH2DT[i] + b32*x2k2[i]);
    }

    for(i=0; i < NPTYPES; i++){
        ptemp[i] = P[i] + h * (b31*DPDT[i] + b32*pk2[i]);
    }    

    

    dynamic(x1temp, x2temp, ptemp, Q1, Q2, CH_H1, CH_H2, CP, x1k3, x2k3, pk3, NH1TYPES, NH2TYPES, NPTYPES);
    
    /*
    for(i=0; i < NTYPES; i++){
        xtemp[i] = X[i] + h * (b41*DXDT[i] + b42*xk2[i] + b43*xk3[i]);
        ytemp[i] = Y[i] + h * (b41*DYDT[i] + b42*yk2[i] + b43*yk3[i]);
    }
    */

    for(i=0; i < NH1TYPES; i++){
        x1temp[i] = XH1[i] + h * (b41*DXH1DT[i] + b42*x1k2[i] + b43*x1k3[i]);
    }

    for(i=0; i < NH2TYPES; i++){
        x2temp[i] = XH2[i] + h * (b41*DXH2DT[i] + b42*x2k2[i] + b43*x2k3[i]);
    }


    for(i=0; i < NPTYPES; i++){
        ptemp[i] = P[i] + h * (b41*DPDT[i] + b42*pk2[i] + b43*pk3[i]);
    }

    dynamic(x1temp, x2temp, ptemp, Q1, Q2, CH_H1, CH_H2, CP, x1k4, x2k4, pk4, NH1TYPES, NH2TYPES, NPTYPES);
    
    /*
    for(i=0; i < NTYPES; i++){
        xtemp[i] = X[i] + h*(b51*DXDT[i] + b52*xk2[i] + b53*xk3[i] + b54*xk4[i]);
        ytemp[i] = Y[i] + h*(b51*DYDT[i] + b52*yk2[i] + b53*yk3[i] + b54*yk4[i]);
    }
    */

    for(i=0; i < NH1TYPES; i++){
        x1temp[i] = XH1[i] + h*(b51*DXH1DT[i] + b52*x1k2[i] + b53*x1k3[i] + b54*x1k4[i]);
    }

    for(i=0; i < NH2TYPES; i++){
        x2temp[i] = XH2[i] + h*(b51*DXH2DT[i] + b52*x2k2[i] + b53*x2k3[i] + b54*x2k4[i]);
    }

    for(i=0; i < NPTYPES; i++){
        ptemp[i] = P[i] + h*(b51*DPDT[i] + b52*pk2[i] + b53*pk3[i] + b54*pk4[i]);
    }

    dynamic(x1temp, x2temp, ptemp, Q1, Q2, CH_H1, CH_H2, CP, x1k5, x2k5, pk5, NH1TYPES, NH2TYPES, NPTYPES);
    
    /*
    for(i=0; i < NTYPES; i++){
        xtemp[i] = X[i] + h*(b61*DXDT[i] + b62*xk2[i] + b63*xk3[i] + b64*xk4[i] + b65*xk5[i]);
        ytemp[i] = Y[i] + h*(b61*DYDT[i] + b62*yk2[i] + b63*yk3[i] + b64*yk4[i] + b65*yk5[i]);
    }
    */

    for(i=0; i < NH1TYPES; i++){
        x1temp[i] = XH1[i] + h*(b61*DXH1DT[i] + b62*x1k2[i] + b63*x1k3[i] + b64*x1k4[i] + b65*x1k5[i]);
    }

    for(i=0; i < NH2TYPES; i++){
        x2temp[i] = XH2[i] + h*(b61*DXH2DT[i] + b62*x2k2[i] + b63*x2k3[i] + b64*x2k4[i] + b65*x2k5[i]);
    }

    for(i=0; i < NPTYPES; i++){
        ptemp[i] = P[i] + h*(b61*DPDT[i] + b62*pk2[i] + b63*pk3[i] + b64*pk4[i] + b65*pk5[i]);
    }

    dynamic(x1temp, x2temp, ptemp,Q1, Q2, CH_H1, CH_H2, CP, x1k6, x2k6, pk6, NH1TYPES, NH2TYPES, NPTYPES);
  
    /*
    for(i=0; i < NTYPES; i++){
        xout[i] = X[i] + h*(c1*DXDT[i] + c3*xk3[i] + c4*xk4[i] + c6*xk6[i]);
        xerr[i] = h*(dc1*DXDT[i] + dc3*xk3[i] + dc4*xk4[i] + dc5*xk5[i] + dc6*xk6[i]);
        yout[i] = Y[i] + h*(c1*DYDT[i] + c3*yk3[i] + c4*yk4[i] + c6*yk6[i]);
        yerr[i] = h*(dc1*DYDT[i] + dc3*yk3[i] + dc4*yk4[i] + dc5*yk5[i] + dc6*yk6[i]);
    }
    */

    for(i=0; i < NH1TYPES; i++){
        x1out[i] = XH1[i] + h * (c1*DXH1DT[i] + c3*x1k3[i] + c4*x1k4[i] + c6*x1k6[i]);
        //std::cout << x1out[i] << std::endl;
        x1err[i] = h*(dc1*DXH1DT[i] + dc3*x1k3[i] + dc4*x1k4[i] + dc5*x1k5[i] + dc6*x1k6[i]);
        //std::cout << x1err[i] << std::endl;
    }

    for(i=0; i < NH2TYPES; i++){
        x2out[i] = XH2[i] + h * (c1*DXH2DT[i] + c3*x2k3[i] + c4*x2k4[i] + c6*x2k6[i]);
        //std::cout << x2out[i] << std::endl;
        x2err[i] = h*(dc1*DXH2DT[i] + dc3*x2k3[i] + dc4*x2k4[i] + dc5*x2k5[i] + dc6*x2k6[i]);
        //std::cout << x2err[i] << std::endl;
    }


    for(i=0; i < NPTYPES; i++){
        pout[i] = P[i] + h * (c1*DPDT[i] + c3*pk3[i] + c4*pk4[i] + c6*pk6[i]);
        perr[i] = h*(dc1*DPDT[i] + dc3*pk3[i] + dc4*pk4[i] + dc5*pk5[i] + dc6*pk6[i]);
    }
}

/***********************************
 * Population dynamics
 ***********************************/
void dynamic(double *XH1, double *XH2, double *P, double **Q1, double **Q2, double *CH_H1, double *CH_H2, double *CP, double *DXH1DT, double *DXH2DT, double *DPDT, int NH1TYPES, int NH2TYPES, int NPTYPES){  
    int i, j;
    //double hostsum[NTYPES], wh[NTYPES];
    double host1sum[NH1TYPES], wh1[NH1TYPES];
    double host2sum[NH2TYPES], wh2[NH2TYPES];
    //double parsum[NTYPES], wp[NTYPES];
    double par1sum[NPTYPES], par2sum[NPTYPES], wp2[NPTYPES];


    //double WH = 0.0;
    double WH1 = 0.0;
    double WH2 = 0.0;
    //double WP = 0.0;
    double WP2 = 0.0;

    double phi1 = PHI;
    double phi2 = (1. - PHI);

    
    /* Calculate interspecies effects */
    /*
    for(i=0; i < NTYPES; i++){
        hostsum[i] = 0;
        for(j=0; j < NTYPES; j++){
            hostsum[i] += Q[i][j]*Y[j];
        }
    }
    */

    /*
    for(j=0; j < NTYPES; j++){
        parsum[j] = 0;
        for(i=0; i < NTYPES; i++){
            parsum[j] += Q[i][j]*X[i];
        }
    }
    */

    for(i=0; i < NH1TYPES; i++){
        host1sum[i] = 0;
        for(j=0; j < NPTYPES; j++){
            host1sum[i] += Q1[i][j]* P[j];
        }
        //std::cout << host1sum[i] << std::endl;
    }

    for(i=0; i < NH2TYPES; i++){
        host2sum[i] = 0;
        for(j=0; j < NPTYPES; j++){
            host2sum[i] += Q2[i][j]* P[j];
        }
    }

    for(j=0; j < NPTYPES; j++){
        par1sum[j] = 0;
        par2sum[j] = 0;
        for(i=0; i < NH1TYPES; i++){
            par1sum[j] += Q1[i][j] * XH1[i];
        }
        for(i=0; i < NH2TYPES; i++){
            par2sum[j] += Q2[i][j] * XH2[i];
        }
    }
     
    /* Fitness functions */
    /*
    for(i=0; i < NTYPES; i++){
        wh[i] = CH[i] * exp(-BETAH*hostsum[i]);
        WH += X[i]*wh[i];
        wp[i] = CP[i] * exp(BETAP*parsum[i]);
        WP += Y[i]*wp[i];
    }
    */

    for(i=0; i < NH1TYPES; i++){
        wh1[i] = CH_H1[i] * exp(-BETAH * host1sum[i]);
        WH1 += XH1[i] * wh1[i];
    }

    for(i=0; i < NH2TYPES; i++){
        wh2[i] = CH_H2[i] * exp(-BETAH * host2sum[i]);
        WH2 += XH2[i] * wh2[i];
    }

    for(i=0; i < NPTYPES; i++){
        wp2[i] = CP[i] *  exp(BETAP * (par1sum[i] * phi1 + par2sum[i] * phi2));
        WP2 += P[i] * wp2[i];
    }
    
    /* Changes in genotype frequencies */
    /*
    for(i=0; i < NTYPES; i++){
        DXDT[i] = ((wh[i]/WH) - 1) * X[i];
    }


    for(i=0; i < NTYPES; i++){
        DYDT[i] = ((wp[i]/WP) - 1) * Y[i];
    }
    */

    for(i=0; i < NH1TYPES; i++){
        DXH1DT[i] = ((wh1[i]/WH1) - 1) * XH1[i];
        //std::cout << "H1 " << DXH1DT[i] << std::endl;
    }

    for(i=0; i < NH2TYPES; i++){
        DXH2DT[i] = ((wh2[i]/WH2) - 1) * XH2[i];
        //std::cout << "H2 " << DXH2DT[i] << std::endl;
    }

    //exit(0);

    for(i=0; i < NPTYPES; i++){
        DPDT[i] = ((wp2[i]/WP2) - 1) * P[i];
        //std::cout << "P " << DPDT[i] << std::endl;
    }
    //exit(0);

}

/***************************************
 * Max and min functions
 ***************************************/
double FMAX(double l, double r)
{
    if( l > r )return l;
    else   return r;
}
double FMIN(double l, double r)
{
    if (l < r )return l;
    else   return r;
}


void initialize_gts(int **GENOS, int NTYPES, int NLOCI, int block){
	
	int j,k,val, switchval,s;
    int switcher = 0;
	
	/* Initialize genotypes */

    for(j=0; j < NLOCI; j++){
        k = 0;
        val = 0;
        switchval = (1 << switcher); // When to switch value between 0 and 1
        while(k < NTYPES){
            s = 0;
            while(s < switchval){
                GENOS[k][j] = val;
                k++;
                s++;
            }
            val = (val+1)%2;
            if(j == block) val = 0;
        }
        if(j != block) switcher ++;
    }
}

/***************************************
 * Function for calculating the number of effective loci
 ***************************************/
void calculateeffective (int NTYPES, int NLOCI, int **GENOTYPES, int *RANGE){
	
	int i, j;
	
	for (i=0; i < NTYPES; i++){
        RANGE[i] = 0;
        for(j=0; j < NLOCI; j++){
            RANGE[i] += GENOTYPES[i][j];
        }
    }
}

/***************************************
 * Function for calculating Q[i][j] and outputting the result
 ***************************************/
 double calculate_Q (int *GENOTYPE_HI, int *GENOTYPE_PJ, double SIGMA, int NLOCI){
	int k;
	int sum;
	std::string hg,pg;
	double Q = 0.0;
	
	hg="";
	pg="";
	sum = 0;
	for (k=0; k < NLOCI; k++){
		hg += std::to_string(GENOTYPE_HI[k]);
		pg += std::to_string(GENOTYPE_PJ[k]);
		if(GENOTYPE_HI[k] > GENOTYPE_PJ[k]) sum += 1;
	}            
	
	Q = pow(SIGMA,sum);
	
	std::cout << "hg" << "\t" << hg << "\tpg:\t" << pg << "\teff\t" << sum << "\tQ=" << Q << std::endl;
	
	return(Q);
}

/***************************************
 * Function to calculate the baseline fitness for each genotype
 ***************************************/
double basefitness (double c_1, double c_2, double normrange){
	
	double bfit = 1.0;
	
	if(c_2 == 0){
		bfit = FMAX(0, 1 - (c_1*normrange) );	
	}else{
		bfit = FMAX(0, 1 - c_1*(1-exp(c_2 * normrange))/(1-exp(c_2)) );
	}
	
	return(bfit);
}


/***************************************
 * Set values for constant arrays
 ***************************************/
void array_values(double *CH_H1, double *CH_H2, double *CP_P2, double **Q1, double **Q2, int *RANGEH1, int *RANGEH2, int *RANGEP, 
                    int NH1TYPES, int NH2TYPES, int NPTYPES, int NLOCI, int fscaleH1, int fscaleH2, int fscaleP, std::string geHfile, std::string geH2file, std::string gePfile){
    
    // double sum;
    int i, j; //, k, s, val, switchval;

    // Allocate memory for an array of pointers
    //int** GENOS = new int*[NTYPES];
    int** GENOS_H1 = new int*[NH1TYPES];
    int** GENOS_H2 = new int*[NH2TYPES];
    int** GENOS_P = new int*[NPTYPES];

    /*
    std::cout << "fscaleP " << fscaleP << std::endl;
    std::cout << "fscaleH1 " << fscaleH1 << std::endl;
    std::cout << "fscaleH2 " << fscaleH2 << std::endl;
    */

    // Allocate memory for each row
    /*
    for (int i = 0; i < NTYPES; ++i) {
        GENOS[i] = new int[NLOCI];
    }
    */

    for(int i = 0; i < NPTYPES; ++i){
        GENOS_P[i] = new int[NLOCI];
    }

    // Allocate memory for each row
    for (int i = 0; i < NH1TYPES; ++i) {
        GENOS_H1[i] = new int[NLOCI];
    }   
 
    // Allocate memory for each row
    for (int i = 0; i < NH2TYPES; ++i) {
        GENOS_H2[i] = new int[NLOCI];
    }   

    //initialize_gts(GENOS, NTYPES, NLOCI, -1);
    //printgenotypes (GENOS, NTYPES, NLOCI);
    initialize_gts(GENOS_P, NPTYPES, NLOCI, blockP);
    //printgenotypes (GENOS_P, NTYPES, NLOCI);


    initialize_gts(GENOS_H1, NH1TYPES, NLOCI, block1);
    //printgenotypes (GENOS_H1, NH1TYPES, NLOCI);
    //std::cout << "Second host" << std::endl;
    initialize_gts(GENOS_H2, NH2TYPES, NLOCI, block2);
    //printgenotypes (GENOS_H2, NH2TYPES, NLOCI);
    //exit(0); 
    
    /* Calculate ranges - use in cost functions */
    //calculateeffective(NTYPES, NLOCI, GENOS, RANGE);
    calculateeffective(NPTYPES, NLOCI, GENOS_P, RANGEP);
    calculateeffective(NH1TYPES, NLOCI, GENOS_H1, RANGEH1);
    calculateeffective(NH2TYPES, NLOCI, GENOS_H2, RANGEH2);

    
    /* Host-parasiteaction matrix, Q */
    /*
    for (i=0; i < NTYPES; i++){
        CH[i] = basefitness(CH1, CH2, (double)RANGE[i]/(double)(fscaleP));	
        CP[i] = basefitness(CP1, CP2, (double)RANGE[i]/(double)(fscaleP));	
        for (j=0; j < NTYPES; j++){
            Q[i][j] = calculate_Q(GENOS[i], GENOS[j], SIGMA, NLOCI);
        }
    }
    */

    for (i=0; i < NPTYPES; i++){
        CP_P2[i] = basefitness(CP1, CP2, (double)RANGEP[i]/(double)(fscaleP));
        // std::cout << "CP_P2[" << i << "]: " <<  CP_P2[i] << std::endl; 	
    }

    /* Host-parasite interaction matrix, Q1 */
    for (i=0; i < NH1TYPES; i++){
        std::cout << (double)RANGEH1[i] << std::endl;
        std::cout << (double)(fscaleH1) << std::endl;
        std::cout << (double)RANGEH1[i]/(double)(fscaleH1) << std::endl;
        CH_H1[i] = basefitness(CH1, CH2, (double)RANGEH1[i]/(double)(fscaleH1));	
        for (j=0; j < NPTYPES; j++){
            Q1[i][j] = calculate_Q(GENOS_H1[i], GENOS_P[j], SIGMA, NLOCI);
        }
    }
    
    std::cout << "Second host" << std::endl;

    /* Host-parasite interaction matrix, Q2 */
    for (i=0; i < NH2TYPES; i++){
        CH_H2[i] = basefitness(CH1, CH2, (double)RANGEH2[i]/(double)(fscaleH2));	
        for (j=0; j < NPTYPES; j++){
            Q2[i][j] = calculate_Q(GENOS_H2[i], GENOS_P[j], SIGMA, NLOCI);
            //std::cout << Q2[i][j] << std::endl;
        }
    }


    //writehostgenotypefile(GENOS, GENOS, GENOS, RANGE, RANGE, CH, CH, Q, Q, NTYPES, NTYPES, NTYPES, NLOCI, geHfile);
    writehostgenotypefile(GENOS_H1, GENOS_H2, GENOS_P, RANGEH1, RANGEH2, CH_H1, CH_H2, Q1, Q2, NH1TYPES, NH2TYPES, NPTYPES, NLOCI, fscaleH1, fscaleH2, geH2file);
    writepathogengenotypefile(GENOS_H1, GENOS_H2, GENOS_P, RANGEP, CP_P2, Q1, Q2, NH1TYPES, NH2TYPES, NPTYPES, fscaleP, gePfile);
    //exit(0);


    for (int i = 0; i < NPTYPES; ++i) {
        //delete[] GENOS[i]; // Free each row
        delete[] GENOS_P[i]; // Free each row
    }
    //delete[] GENOS; // Free the array of pointers

    for (int i = 0; i < NH1TYPES; ++i) {
        delete[] GENOS_H1[i]; // Free each row
    }

    for (int i = 0; i < NH2TYPES; ++i) {
        delete[] GENOS_H2[i]; // Free each row
    }

    delete[] GENOS_P; // Free the array of pointers
    delete[] GENOS_H1; // Free the array of pointers
    delete[] GENOS_H2; // Free the array of pointers


}


void extinct_check(double *FREQS, int NTYPES, double minfreq){
    for(int i=0; i < NTYPES; i++){
        FREQS[i] = (FREQS[i] > 1e-6) ? FREQS[i] : 0.0;
    }
    //std::cout << "Done with extinction check" << std::endl;
}

