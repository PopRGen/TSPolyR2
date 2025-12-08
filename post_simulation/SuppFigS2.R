suppressPackageStartupMessages(require("tidyverse"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("patchwork"))
suppressPackageStartupMessages(require("ggpubr"))




plotdir <- "Figures/supplementary_figures"


if(!dir.exists(plotdir)){
  dir.create(plotdir)
}


# These parameters are fixed throughout the simulations
sigma <- 0.85
betaH <- 1
betaP <- 1

# Settings for the simulation plotted in Fig1 D and E 
prefix <- "ext_reps"
suffix <- "2025-01-07"
omegaH <- 0.1
omegaP <- 0.1
XiH <- 3
XiP <- 3
phi <- 0.5
seed <- 1600
plotdir <- "Figures/main_figures"
indir <- "Data_SFig2"


if(!dir.exists(plotdir)){
  dir.create(plotdir, recursive = TRUE)
}

# These parameters are fixed throughout the simulations
sigma <- 0.85
betaH <- 1
betaP <- 1


# Simulation prefix. This is the entire file name prefix for all the result files of the simulation
simprefix <- paste0(prefix,
                    "_phi_",phi,
                    "_cH_1_",omegaH,
                    "_cH_2_",XiH,
                    "_cP_1_",omegaP,
                    "_cP_2_",XiP,
                    "_s_",sigma,
                    "_bH_",betaH,
                    "_bP_",betaP,
                    "_seed_",seed,
                    "_",suffix)

# Simulation prefix. This is the entire file name prefix for all the result files of the simulation
simprefix <- paste0(indir,"/",prefix,
                    "_phi_", phi,
                    "_cH_1_", omegaH,
                    "_cH_2_", XiH,
                    "_cP_1_", omegaP,
                    "_cP_2_", XiP,
                    "_s_", sigma,
                    "_bH_", betaH,
                    "_bP_", betaP,
                    "_seed_", seed,
                    "_", suffix)

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
  ylim(0,1) + 
  theme_bw(base_size = 14)  + 
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
  theme_bw(base_size = 14)  + 
  xlab(expression(P[110])) + 
  ylab(expression(M["001"])) + 
  ggtitle("Host M")


phase_3 <- ggplot(daten_sub , aes(P_3, H1_1)) + 
  geom_path(color = "gray50") +
  geom_point() +
  geom_point(data = daten_sub |> head(n=1), aes(P_3, H1_1), size = 3, shape = 17) +
  geom_point(data = daten_sub |> tail(n=1), aes(P_3, H1_1), size = 3, shape = 15) +
  xlim(0,0.5) +
  ylim(0,1) + 
  theme_bw(base_size = 14)  + 
  xlab(expression(P[110])) + 
  ylab(bquote(H["010"]))


phase_4 <- ggplot(daten_sub , aes(P_3, H2_1)) + 
  geom_path(color = "gray50") +
  geom_point() +
  geom_point(data = daten_sub |> head(n=1), aes(P_3, H2_1), size = 3, shape = 17) +
  geom_point(data = daten_sub |> tail(n=1), aes(P_3, H2_1), size = 3, shape = 15) +
  xlim(0,0.5) +
  ylim(0,1) + 
  theme_bw(base_size = 14)  + 
  xlab(expression(P[110])) + 
  ylab(expression(M["100"]))


SFig2 <- (phase_1 + phase_2) / (phase_3 + phase_4) + 
  plot_layout(guides = "collect") & 
  # Add labels to the subfigures
  plot_annotation(
    tag_levels = list(c("(a)", "(b)", "(c)", "(d)")),
    tag_prefix = "",
    tag_suffix = ""
  ) &
  # Define the size of the subfigure labels and the plot subplot titles
  theme(
    plot.tag = element_text(face = 2, size = 14),
    plot.title = element_text(face = 2, size = 16)
  )

pdf(paste0(plotdir,"/SFig2_9120_9250_pre.pdf"), width = 6.5, height = 6)
print(SFig2)
dev.off()
