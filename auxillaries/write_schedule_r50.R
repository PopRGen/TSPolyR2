require("tidyverse")
require("argparser")

# Script to generate a simulation schedule

# Set up command line argument parser
args <- arg_parser("Parser to generate simulation schedule")
args <- add_argument(args,"--executable", 
                     help="Path to the executable to use", 
                     default = "Non_ecological_clean")
args <- add_argument(args, "--prefix",
                     help="Prefix to append to the output files", 
                     default = "ext_reps")
args <- add_argument(args, "--suffix", 
                     help="suffix to append to the output files", 
                     default = format(Sys.Date(),"%Y-%m-%d"))
args <- add_argument(args, "--cH_2", 
                     help="cH_2 value to use", 
                     default = 3,
                     type = "numeric")
args <- add_argument(args, "--cP_2", 
                     help="cP_2 value to use", 
                     default = 0,
                     type = "numeric")
args <- add_argument(args, "--phi", 
                     help="phi value for the simulations", 
                     default = 0.5,
                     type = "numeric")
args <- add_argument(args, "--nrep", 
                     help="number of different initial seed", 
                     default = 50,
                     type = "numeric")
args <- add_argument(args, "--seed_start", 
                     help="lower seed to start with", 
                     default = 400,
                     type = "numeric")
args <- add_argument(args, "--seed_interval", 
                     help="seed interval tp ise", 
                     default = 400,
                     type = "numeric")
args <- add_argument(args, "--outdir", 
                     help="directory to which the sheet should be placed", 
                     default = "schedules_r50",
                     type = "numeric")


# Parse the command line arguments
pargs <- parse_args(args)
executable <- pargs[["executable"]]
prefix <- pargs[["prefix"]]
suffix <- pargs[["suffix"]]
CH2 <- pargs[["cH_2"]]
CP2 <- pargs[["cP_2"]]
phi <- pargs[["phi"]]
nrep <- pargs[["nrep"]]
outdir <- pargs[["outdir"]]
seed_start <- pargs[["seed_start"]]
seed_interval <- pargs[["seed_interval"]]

outrandom <- paste0(outdir,"/schedule_",prefix,"_CH2_",CH2,"_CP2_",CP2,"_phi_",phi,"_",suffix,"_random.txt")

# Set up the small fixed schedule
seed <- seq(seed_start, seed_start+ seed_interval*(nrep-1), by=seed_interval)

# Write simulation schedule for simulations with random initial frequencies
daten_random <- expand.grid(phi=phi, cH_2=CH2, cP_2=CP2, seed=seed)

# Write the simulation schedule for the simulations with random initial frequencies
schedule_random <- daten_random %>%
  mutate(executable=executable,
         prefix=prefix,
         suffix=suffix,
         init="r") %>%
  relocate(executable)

# Write the result to file
write_delim(schedule_random, outrandom,col_names = FALSE)
