require("tidyverse")
require("argparser")


args <- arg_parser("Command line argument parser")
args <- add_argument(args,arg="--prefix",
                     help="Prefix to append to all simulations",
                     default = "Non_ecological",
                     type="character")
args <- add_argument(args,arg="--suffix",
                     help="Suffix to append to all filenames",
                     type="character",
                     default= format(Sys.Date(),"%Y-%m-%d"))
args <- add_argument(args,arg="--executable",
                     help="Name of the compiled C++-code used for running the simulations",
                     default ="Non_ecological_clean",
                     type="character")
args <- add_argument(args,arg="--cH_2",
                     help="Value of cH_2 to use for the simulations",
                     default ="3",
                     type="character")
args <- add_argument(args,arg="--cP_2",
                     help="Value of cP_2 to use for simulations",
                     default ="3",
                     type="character")
args <- add_argument(args, arg = "--block1",
                     help="Which locus to block of in host1. Note is zero-based",
                     default=0,
                     type = "numeric")
args <- add_argument(args, arg = "--block2",
                     help="Which locus to block of in host2. Note is zero-based",
                     default=1,
                     type = "numeric")
args <- add_argument(args,arg="--sigma",
                     help="Value of sigma to use for simulations",
                     default = 0.85,
                     type="numeric")
args <- add_argument(args,arg="--betaH",
                     help="Value of betaH for simulations",
                     default = 1,
                     type="numeric")
args <- add_argument(args,arg="--betaP",
                     help="Value of betaP to use for simulations",
                     default = 1,
                     type="numeric")
args <- add_argument(args,arg="--nloci",
                     help="Number of loci in the simulation",
                     default = 3,
                     type="numeric")
args <- add_argument(args,arg="--init",
                     help="Initial condition",
                     default = "f",
                     type="character")
args <- add_argument(args,arg="--excheck",
                     help="Check if genotypes went extinct. 0= no, 1 = yes",
                     default = 0,
                     type="numeric")
args <- add_argument(args,arg="--cH1_min",
                     help="Lower limit to draw the CH1 values from",
                     default = 0.1,
                     type="numeric")
args <- add_argument(args,arg="--cP1_min",
                     help="Lower limit to draw the CP1 values from",
                     default = 0.1,
                     type="numeric")
args <- add_argument(args,arg="--cH1_max",
                     help="Upper limit for the the uniform from which the cH1 values are to be drawn",
                     default = 0.3,
                     type="numeric")
args <- add_argument(args,arg="--cP1_max",
                     help="Upper limit for the the uniform from which the cP1 values are to be drawn",
                     default = 0.3,
                     type="numeric")
args <- add_argument(args,arg="--phi_min",
                     help="Upper limit for the the uniform from which the cH1 values are to be drawn",
                     default = 0.5,
                     type="numeric")
args <- add_argument(args,arg="--phi_max",
                     help="Upper limit for the the uniform from which the cP1 values are to be drawn",
                     default = 0.5,
                     type="numeric")
args <- add_argument(args,arg="--rseed",
                     help="Seed for the RNG of the Rscript",
                     default = 2,
                     type="numeric")
args <- add_argument(args,arg="--outdir",
                     help="Output directory for the simulation results",
                     default = "../../results_random",
                     type="character")
args <- add_argument(args,arg="--nreps",
                     help="Output directory for the simulation results",
                     default = 10,
                     type="numeric")


####################################
# MAIN ######################
##############################

# Parse command line arguments
pargs <- parse_args(args)

prefix <- pargs[["prefix"]]
suffix <- pargs[["suffix"]]
executable <- pargs[["executable"]]
dir_exec <- dirname(executable)
cH_2 <- as.numeric(pargs[["cH_2"]])
cP_2 <- as.numeric(pargs[["cP_2"]])
phi <- as.numeric(pargs[["phi"]])
block1 <- as.numeric(pargs[["block1"]])
block2 <- as.numeric(pargs[["block2"]])
sigma <- as.numeric(pargs[["sigma"]])
betaH <- as.numeric(pargs[["betaH"]])
betaP <- as.numeric(pargs[["betaP"]])
nloci <- as.numeric(pargs[["nloci"]])
rseed <- as.numeric(pargs[["rseed"]])
excheck <- as.numeric(pargs[["excheck"]])
init <- pargs[["init"]]
outdir <- pargs[["outdir"]]
cH1_min <- pargs[["cH1_min"]]
cP1_min <- pargs[["cP1_min"]]
cP2_min <- pargs[["cH2_min"]]
cP2_min <- pargs[["cP2_min"]]
cH1_max <- pargs[["cH1_max"]]
cP1_max <- pargs[["cP1_max"]]
cP2_max <- pargs[["cH2_max"]]
cP2_max <- pargs[["cP2_max"]]
phi_min <- pargs[["phi_min"]]
phi_max <- pargs[["phi_max"]]
nreps <- pargs[["nreps"]]

# Set the seed for the simulations
set.seed(rseed)

# Assign to CH2 and CP2 for historic reasons
CH2 <- cH_2
CP2 <- cP_2


# Generate the output directory structure
if(!dir.exists(outdir)){
  dir.create(outdir)
}


outdir_simulations <- paste0(outdir,"/",prefix,"_CH2_",CH2,"_CP2_",CP2,"_random_combis_", suffix,"/simulations")
#outdir_simulations <- "/scratch/hm2840/modelling_rerun/results_random/simulations_CH2_3_CP2_3_random_combis"

if(!dir.exists(outdir_simulations)){
  dir.create(outdir_simulations, recursive = TRUE)
}




combis <- data.frame(cH_1 = round(runif(nreps, cH1_min, cH1_max), digits = 5),
                     cP_1 = round(runif(nreps, cP1_min, cP1_max), digits = 5),
                     phi = round(runif(nreps, phi_min, phi_max), digits = 5),
                     seed = sample(1:25000, nreps))

##################################
## MAIN ##########################
##################################

fprefix <- paste0("CH2_", CH2, "_CP2_", CP2,"_random_r_",nreps,"_RSEED_",rseed)

write_tsv(combis,  paste0(outdir_simulations,"/", fprefix,".tsv"))

tararchive <- paste0(outdir_simulations, "/", fprefix, ".tar")
system(paste0("cd ",outdir_simulations,";", paste("tar -cf", tararchive, paste0(fprefix, ".tsv"))))

outfile <- paste0(outdir_simulations, "/", fprefix, "_end.txt")
if(file.exists(outfile)){system(paste("rm", outfile))}
system(paste("touch", outfile))

for(i in 1:nrow(combis)){
  
  CH1 <- round(combis[i,"cH_1"], digits = 3)
  CP1 <- round(combis[i,"cP_1"], digits = 3)
  phi <- round(combis[i,"phi"], digits = 3)
  seed <- combis[i,"seed"]
  simprefix <- paste0(outdir_simulations,"/", fprefix,"_simulation_", sprintf("%05d", i))
    
  print(paste("Before running simulation:", sprintf("%05d",i)))
    
  # Run the actual simulation
  system(paste(executable,
                 "--phi",phi,
                 "--cH_1",CH1,
                 "--cH_2",CH2,
                 "--cP_1",CP1,
                 "--cP_2",CP2,
                 "--sigma",sigma,
                 "--betaH",betaH,
                 "--betaP",betaP,
                 "--block1",block1,
                 "--block2", block2,
                 "--nloci",nloci,
                 "--init", init,
                 "--fprefix",simprefix,
                 "--excheck",excheck,
                 "--seed",seed))
  print(paste("After running simulation:", sprintf("%05d",i)))
    
  last_line <- system(paste("tail -n1",  paste0(outdir_simulations,"/", fprefix,"_simulation_", sprintf("%05d", i), "_dynamics.txt")), intern = T)
  system(paste0("echo ",last_line,"simulation_", sprintf("%05d", i), " >>", outfile))
    
  if( (i%%50) == 0 ){
    system(paste0("cd ",outdir_simulations,";", paste("tar -rf", tararchive, paste0(fprefix, "_simulation*"))))
    system(paste0("cd ",outdir_simulations,"; rm *_simulation*"))
  }
}
