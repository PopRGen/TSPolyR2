suppressPackageStartupMessages(require("tidyverse"))
suppressPackageStartupMessages(require("argparser"))
suppressPackageStartupMessages(require("paletteer"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("patchwork"))
suppressPackageStartupMessages(require("ggpubr"))


args <- arg_parser("Command line argument parser")
args <- add_argument(args,arg="--prefix",
                     help="Prefix to append to all simulations",
                     default = "ext_reps",
                     type="character")
args <- add_argument(args,arg="--suffix",
                     help="Suffix to append to all filenames",
                     type="character",
                     default="2025-01-07")
args <- add_argument(args,arg="--phi",
                     help="Value of phi to use for the simulations",
                     default =0.5,
                     type="numeric")
args <- add_argument(args,arg="--cH_2",
                     help="Value of cH_2 to use for the simulations",
                     default =3,
                     type="numeric")
args <- add_argument(args,arg="--cP_2",
                     help="Value of cP_2 to use for simulations",
                     default =3,
                     type="numeric")
args <- add_argument(args,arg="--cH_1",
                     help="Value of cH_1 to use for the simulations",
                     default = 0.1,
                     type="numeric")
args <- add_argument(args,arg="--cP_1",
                     help="Value of cP_1 to use for simulations",
                     default = 0.1,
                     type="numeric")
args <- add_argument(args,arg="--parent_outdir",
                     help="Top level directory for results",
                     default ="results_r50",
                     type="character")
args <- add_argument(args,arg="--seed",
                     help="Seed to use for the simulation",
                     default = 1600,
                     type="numeric")
args <- add_argument(args, arg="--plotdir",
                     help="Directory for storing the simulations",
                     default="Figures_seed_comparison",
                     type="character")

# Parse command line arguments
pargs <- parse_args(args)

prefix <- pargs[["prefix"]]
suffix <- pargs[["suffix"]]
CH1 <- as.numeric(pargs[["cH_1"]])
CP1 <- as.numeric(pargs[["cP_1"]])
CH2 <- as.numeric(pargs[["cH_2"]])
CP2 <- as.numeric(pargs[["cP_2"]])
phi <- as.numeric(pargs[["phi"]])
seed <- as.numeric(pargs[["seed"]])
parent_outdir <- pargs[["parent_outdir"]]
plotdir <- pargs[["plotdir"]]
# 

plotdir <- paste0(plotdir)


if(!dir.exists(plotdir)){
  dir.create(plotdir)
}


# These parameters are fixed throughout the simulations
sigma <- 0.85
betaH <- 1
betaP <- 1


# Find the directory with the corresponding CH2 and CP2 combination
outdir_cHcP <- paste0(parent_outdir,"/",prefix,"_CH2_",CH2,"_CP2_",CP2,"_",suffix)
outdir_simulations <- paste0(outdir_cHcP,"/simulations/simulations_CH2_",CH2,"_CP2_",CP2,"_",seed)



# Simulation prefix. This is the entire file name prefix for all the result files of the simulation
simprefix <- paste0(outdir_simulations,"/",prefix,"_phi_",phi,"_cH_1_",CH1,"_cH_2_",CH2,"_cP_1_",CP1,"_cP_2_",CP2,"_s_",sigma,"_bH_",betaH,"_bP_",betaP,"_seed_",seed,"_",suffix)

# Read the dynamics from the file
daten <- read.table(paste0(simprefix,"_dynamics.txt"))
names(daten) <- c("time",paste0("H1_",0:3), paste0("H2_",0:3), paste0("P_",0:7))

# Convert the dynamics to a tibble
daten <- as_tibble(daten) 

# Get the information on the genotypes for the hosts
hl_info <- read.table(paste0(simprefix,"_genotypeH2.txt"),header=T,
                      colClasses=c("character","character","character",rep("numeric",times=10))) %>%
  as_tibble() %>%
  mutate(HostID=paste(Host,ID,sep="_"))

# Read the pathogen genotype info file
pl_info <- read.table(paste0(simprefix,"_genotypeP.txt"),header=T,
                      colClasses=c("character","character","character", rep("numeric",times=10))) %>%
  as_tibble() %>%
  mutate(ParaID=paste(Para,ID,sep="_"))


daten_sub <- daten |> filter(time > 9120 & time < 9250)

phase_1 <- ggplot(daten_sub , aes(P_3, H1_2)) + 
  geom_path(color = "gray50") +
  geom_point() +
  geom_point(data = daten_sub |> head(n=1), aes(P_3, H1_2), size = 3, shape = 17) +#+
  geom_point(data = daten_sub |> tail(n=1), aes(P_3, H1_2), size = 3, shape = 15) +
  xlim(0,0.5) +
  ylim(0,1) + theme_bw()  + 
  xlab(expression(P[110])) + 
  ylab(expression(H["001"])) +
  ggtitle("Host H")


phase_2 <- ggplot(daten_sub , aes(P_3, H2_2)) + 
  geom_path(color = "gray50") +
  geom_point() +
  geom_point(data = daten_sub |> head(n=1), aes(P_3, H2_2), size = 3, shape = 17) +
  geom_point(data = daten_sub |> tail(n=1), aes(P_3, H2_2), size = 3, shape = 15) +
  xlim(0,0.5) +
  ylim(0,1) + 
  theme_bw() + 
  xlab(expression(P[110])) + 
  ylab(expression(M["001"])) + 
  ggtitle("Host M")


phase_3 <- ggplot(daten_sub , aes(P_3, H1_1)) + 
  geom_path(color = "gray50") +
  geom_point() +
  geom_point(data = daten_sub |> head(n=1), aes(P_3, H1_1), size = 3, shape = 17) +
  geom_point(data = daten_sub |> tail(n=1), aes(P_3, H1_1), size = 3, shape = 15) +
  xlim(0,0.5) +
  ylim(0,1) + theme_bw()  + 
  xlab(expression(P[110])) + 
  ylab(bquote(H["010"]))


phase_4 <- ggplot(daten_sub , aes(P_3, H2_1)) + 
  geom_path(color = "gray50") +
  geom_point() +
  geom_point(data = daten_sub |> head(n=1), aes(P_3, H2_1), size = 3, shape = 17) +
  geom_point(data = daten_sub |> tail(n=1), aes(P_3, H2_1), size = 3, shape = 15) +
  xlim(0,0.5) +
  ylim(0,1) + theme_bw() + 
  xlab(expression(P[110])) + 
  ylab(expression(M["100"]))


suppl1 <- (phase_1 + phase_2) / (phase_3 + phase_4) + 
  plot_layout(guides = "collect")

pdf(paste0(plotdir,"/supplement_Fig1_9120_9250.pdf"), width = 6.5, height = 6)
print(suppl1)
dev.off()
