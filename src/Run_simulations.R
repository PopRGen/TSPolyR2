require("tidyverse")
library("openxlsx")
require("xtable")
library("argparser")


source("Run_simulations_functions.R")

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
args <- add_argument(args,arg="--phi",
                     help="Value of phi to use for simulations",
                     default ="0.5",
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
args <- add_argument(args,arg="--parent_outdir",
                     help="Top level directory for results",
                     default ="results",
                     type="character")
args <- add_argument(args,arg="--seed",
                     help="Seed to use for the simulation",
                     default ="1",
                     type="character")
args <- add_argument(args,arg="--nloci",
                     help="Number of loci in the simulation",
                     default = 3,
                     type="numeric")
args <- add_argument(args,arg="--init",
                     help="Initial condition",
                     default = "f",
                     type="character")
args <- add_argument(args,arg="--excheck",
                     help="Initial condition",
                     default = 0,
                     type="numeric")
args <- add_argument(args,arg="--gsize",
                     help="For which grid should simulations be run. standard: cH_1=0.1 , cP_1=0.1,
                     small: cH_1={0.1,0.2,0.3} and cP_1={0.1,0.2,0.3}, full: cH_1={seq(0.05,0.35,by=0.05)} and cP_1={seq(0.05,0.35,by=0.05)}, 
                     additional: complements small to obtain full grid",
                     default ="full",
                     type="character")




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
seed <- as.numeric(pargs[["seed"]])
excheck <- as.numeric(pargs[["excheck"]])
init <- pargs[["init"]]
parent_outdir <- pargs[["parent_outdir"]]
gsize <- pargs[["gsize"]]

if(!gsize%in%c("standard","small","full", "additional")){stop("--gsize must be one of standard, small, full")}

# Parameter settings

if(gsize%in%c("full", "additional")){
  cH_1 <- seq(0.05,0.35,len=7)
  cP_1 <- seq(0.05,0.35,len=7)
}else{
  if(gsize=="small"){
    cH_1 <- seq(0.1,0.3,len=3)
    cP_1 <- seq(0.1,0.3,len=3)
  }else{
    cH_1 <- seq(0.1)
    cP_1 <- seq(0.1)
  }
}


# Assign to CH2 and CP2 for historic reasons
CH2 <- cH_2
CP2 <- cP_2


# Generate the output directory structure
if(!dir.exists(parent_outdir)){
  dir.create(parent_outdir)
}

outdir_cHcP <- paste0(parent_outdir,"/",prefix,"_CH2_",CH2,"_CP2_",CP2,"_",suffix)
outdir_simulations <- paste0(outdir_cHcP,"/simulations/simulations_CH2_",CH2,"_CP2_",CP2,"_",seed)
outdir <- outdir_cHcP

if(!dir.exists(outdir_cHcP)){
  dir.create(outdir_cHcP,recursive = TRUE)
}

if(!dir.exists(outdir_simulations)){
  dir.create(outdir_simulations,recursive = TRUE)
}


baseoutname <- paste0(prefix,"_phi_",phi,"_cH_2_",CH2,"_cP_2_",CP2,"_s_",sigma,"_bH_",betaH,"_bP_",betaP,"_seed_",seed)


##################################
## MAIN ##########################
##################################


for(CH1 in cH_1){
  for(CP1 in cP_1){
    if((((CH1*20)%%2 + (CP1*20)%%2) == 0) & gsize == "additional"){break}
    # Simulation prefix
    simprefix <- paste0(outdir_simulations,"/",
                        prefix,"_phi_",phi,
                        "_cH_1_",CH1,
                        "_cH_2_",CH2,
                        "_cP_1_",CP1,
                        "_cP_2_",CP2,
                        "_s_",sigma,
                        "_bH_",betaH,
                        "_bP_",betaP,
                        "_seed_",seed,
                        "_",suffix)
    # Prepare the labels
    plabel <- paste0(prefix,"_phi_",phi,
                     "_cH_1_",CH1,
                     "_cH_2_",CH2,
                     "_cP_1_",CP1,
                     "_cP_2_",CP2,
                     "_s_",sigma,
                     "_bH_",betaH,
                     "_bP_",betaP,
                     "_seed_",seed,
                     "_",suffix)
    hlabel <- paste0("phi=",phi,"     sigma=",sigma,
                     "\ncH_1=",CH1,
                     "  cP_1=",CP1,
                     "\ncH_2=",CH2,"     cP_2=",
                     CP2,"\nbetaH=",betaH,
                     "    betaP=",betaP,
                     "\nseed=",seed)
    
    
    print("Before running simulation")
    
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
    print("After running simulation")
    
  }
}
