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
                     default="Figures",
                     type="character")

plotrange <- function(daten_longer_cumsum, tmin, tmax){
  HRange <- ggplot(daten_longer_cumsum |> filter(time  > tmin & time < tmax & freq > 0 & Species%in%c("H1", "H2"))) +
    geom_segment(aes(x = time, xend = time, y = rsortID, yend = rsortID + 1, color = Genotype, alpha = freq)) +
    scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                   "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20"), drop= F) +
    scale_alpha_identity() +
    facet_grid(rows = vars(Species)) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    scale_y_continuous(breaks=c(0.5, 2, 3.5), labels=c("R0","R1","R2"), limits = c(0,4), expan = expansion()) +
    scale_x_continuous(limits = c(tmin, tmax), expan = expansion()) +
    geom_hline(yintercept = c(1, 3), color = 'gray') +
    ylab("Range") + 
    xlab("")
  
  PRange <-  ggplot(daten_longer_cumsum |> filter(time  > tmin & time < tmax & freq > 0 & Species%in%c("P"))) +
    geom_segment(aes(x = time, xend = time, y = rsortID, yend = rsortID + 1, color = Genotype, alpha = freq)) +
    scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                   "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
    scale_alpha_identity() +
    facet_grid(rows = vars(Species)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_continuous(breaks=c(0.5, 2.5, 5.5, 7.5), labels=c("R0","R1","R2","R3"), 
                       limits = c(0,8), expan = expansion()) +
    scale_x_continuous(limits = c(tmin, tmax), expan = expansion()) +
    geom_hline(yintercept = c(1, 4, 7), color = 'gray') +
    ylab("Range")
  
  Rangeplt <- ggarrange(HRange, PRange, legend.grob = get_legend(PRange), nrow = 2, legend = "right")
  return(Rangeplt)
}


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

outplt <- paste0(plotdir,"/dynamics_cH2_",CH2,"_cP2_",CP2,"_phi_",phi,"_cH1_",CH1,"_cP1_",CP1,".pdf")
outplt_extended <- paste0(plotdir,"/dynamics_cH2_",CH2,"_cP2_",CP2,"_phi_",phi,"_cH1_",CH1,"_cP1_",CP1,"_extended.pdf")

if(!dir.exists(plotdir)){
  dir.create(plotdir, recursive = TRUE)
}


# These parameters are fixed throughout the simulations
sigma <- 0.85
betaH <- 1
betaP <- 1


# Find the directory with the corresponding CH2 and CP2 combination
outdir_cHcP <- paste0(parent_outdir,"/",prefix,"_CH2_",CH2,"_CP2_",CP2,"_",suffix)
outdir_simulations <- paste0(outdir_cHcP,"/all_dynamics_long")



# Simulation prefix. This is the entire file name prefix for all the result files of the simulation
simprefix <- paste0(prefix,"_phi_",phi,"_cH_1_",CH1,"_cH_2_",CH2,"_cP_1_",CP1,"_cP_2_",CP2,"_s_",sigma,"_bH_",betaH,"_bP_",betaP,"_seed_",seed,"_",suffix)

# Read the dynamics from the file
daten_longer_cumsum <- read_tsv(paste0(simprefix,"_long.tsv"))

daten_longer_cumsum <- daten_longer_cumsum |>
  ungroup() |>
  mutate(Species_label = case_when(
    Species == "H1" ~ "Host H",
    Species == "H2" ~ "Host M",
    Species == "P" ~ "Pathogen",
    TRUE ~ "unidentified"
  ))

# Add additional column to allow for proper formatting of the final plot
daten_longer_cumsum <- daten_longer_cumsum |>
  mutate(IPart = case_when(
    Species%in%c("H1", "H2") ~ "Hosts",
    Species == "P" ~ "Pathogen",
    TRUE ~ "unidentified"
  )) |>
  mutate(Genotype = factor(Genotype, levels = c("000", "100", "010", "001", "110", "101", "011", "111")))


dynamicflt_plt_full <- ggplot(daten_longer_cumsum,
                              aes(x=time)) +
  annotate("rect", fill = "gray20", alpha = 0.5, xmin = c(0, 9000, 99000), xmax = c(1000, 10000, 100000), ymin = -Inf, ymax = Inf) +  
  geom_line(aes(x=time, y=freq, color = Genotype),  linewidth = 0.2) +
  scale_alpha_manual(values=c(1, 0.5, 1)) +
  facet_grid(rows = vars(Species_label), scales = "free_y") +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.x = element_line(color = "gray40"), panel.grid.minor.x = element_line(color = "gray80")) +
  ylab("freq") +
  #labs(subtitle = "(a)") +
  #scale_linetype_manual(values = c("Host H"="solid", "Host M"="11", "Pathogen"="solid")) +
  ylab("Genotype frequency") +
  ylim(0,1)  +
  theme(plot.margin = unit(c(0.05, 0.2, 0, 1), 'lines')) +
  scale_x_continuous(limits = c(0, 100000), expan = expansion())   +
  theme(legend.key.width = unit(1,"cm"),
        #legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(hjust = 0.5, size = 12, face = 2),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        panel.grid.major.x = element_line(color = "gray80"),
        axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm"), size = 16), #+#,
        axis.text.x = element_text(size = 14)) + #,
        #axis.title = element_text(size = 12))+ 
  guides(color = "none", linetype = "none", alpha = "none", shape = "none") 
  
svg("fig1_full.svg", width = 14, height = 4)
print(dynamicflt_plt_full)
dev.off()


dynamicflt_plt1 <- ggplot(daten_longer_cumsum |>
                           filter(time < 1000 & freq > 0),
                         aes(x=time)) +
  geom_line(aes(x=time, y=freq, linetype = Species_label, color = Genotype, alpha = Species_label)) +
  scale_alpha_manual(values=c(1, 0.5, 1)) +
  #geom_line(aes(x=time, y=freq, color = Genotype, group = interaction(Species, Genotype))) +
  geom_point(aes(x=time, y=freq, shape = Species_label, color = Genotype, alpha = Species_label), size = 0.4) +
  #geom_line(aes(y = cumsum_freq, color = Genotype)) +
  facet_grid(rows = vars(IPart), scales = "free_y") +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  theme_bw() +
  theme(panel.grid.major.x = element_line(color = "gray40"), panel.grid.minor.x = element_line(color = "gray80")) +
  ylab("freq") +
  #labs(subtitle = "(a)") +
  scale_linetype_manual(values = c("Host H"="solid", "Host M"="11", "Pathogen"="solid")) +
  ylab("Frequency") +
  ggtitle("time span: 0 - 1,000") +
  ylim(0,1)  +
  theme(plot.margin = unit(c(0.05, 0.2, 0, 1), 'lines')) +
  scale_x_continuous(limits = c(0, 1000), expan = expansion())   +
  theme(legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(hjust = 0, size = 12, face = 2),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))+ 
  guides(color = "none", linetype = "none", alpha = "none", shape = "none")


dynamicflt_plt2 <- ggplot(daten_longer_cumsum |>
         filter(time < 10000 & time > 9000 & freq > 0),
       aes(x=time)) +
  geom_line(aes(x=time, y=freq, linetype = Species_label, color = Genotype, alpha = Species_label)) +
  scale_alpha_manual(values=c(1, 0.5, 1)) +
  #geom_line(aes(x=time, y=freq, color = Genotype, group = interaction(Species, Genotype))) +
  geom_point(aes(x=time, y=freq, shape = Species_label, color = Genotype, alpha = Species_label), size = 0.4) +
  #geom_line(aes(y = cumsum_freq, color = Genotype)) +
  facet_grid(rows = vars(IPart), scales = "free_y") +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  theme_bw() +
  theme(panel.grid.major.x = element_line(color = "gray20"), panel.grid.minor.x = element_line(color = "gray80")) +
  ylab("freq") +
  #labs(subtitle = "(b)") +
  scale_linetype_manual(values = c("Host H"="solid", "Host M"="11", "Pathogen"="solid")) +
  ylab("Frequency") +
  ggtitle("time span: 9,000 - 10,000") +
  ylim(0,1) +
  theme(plot.margin = unit(c(0.05,0.2,0,1), 'lines')) + 
  scale_x_continuous(limits = c(9000, 10000), expan = expansion())    +
  theme(legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(hjust = 0, size = 12, face = 2),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) +
  guides(color = "none", alpha = "none", 
         linetype = guide_legend(override.aes = list(color = "gray50", size = c(4,4,4), linewidth = c(2,1,2), alpha = c(1,0.8,1))))


dynamicflt_plt3 <- ggplot(daten_longer_cumsum |>
                           filter(time > 99000 & time < 100000 & freq > 0),
                         aes(x=time)) +
  geom_line(aes(x=time, y=freq, linetype = Species_label, color = Genotype, alpha = Species_label)) +
  scale_alpha_manual(values=c(1, 0.5, 1)) +
  #geom_line(aes(x=time, y=freq, color = Genotype, group = interaction(Species, Genotype))) +
  geom_point(aes(x=time, y=freq, shape = Species_label, color = Genotype, alpha = Species_label), size = 0.4) +
  #geom_line(aes(y = cumsum_freq, color = Genotype)) +
  facet_grid(rows = vars(IPart), scales = "free_y") +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  theme_bw() +
  theme(panel.grid.major.x = element_line(color = "gray40"), panel.grid.minor.x = element_line(color = "gray80")) +
  ylab("freq") + #+
  #labs(subtitle = "(c)") +
  scale_linetype_manual(values = c("Host H"="solid", "Host M"="11", "Pathogen"="solid")) +
  ylab("Frequency") +
  ylim(0,1) +
  ggtitle("time span: 99,000 - 100,000") +
  theme(plot.margin = unit(c(0.05, 0.2,0,1), 'lines')) +
  scale_x_continuous(limits = c(99000, 100000), expan = expansion()) +
  theme(legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(hjust = 0, size = 12, face = 2),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) + 
  guides(color = "none", linetype = "none", alpha = "none", shape = "none")


daten_longer_cumsum <- daten_longer_cumsum |>
  separate_wider_position(Genotype, widths = c("L1"=1, "L2"=1, "L3"=1) , cols_remove = F) |>
  mutate(range = as.numeric(L1) + as.numeric(L2) + as.numeric(L3)) |>
  mutate(rsortID = case_when(
    Species %in% c("H1","H2") ~ as.numeric(ID),
    Species == "P" & Genotype == "000" ~ 0,
    Species == "P" & Genotype == "100" ~ 1,
    Species == "P" & Genotype == "010" ~ 2,
    Species == "P" & Genotype == "001" ~ 3,
    Species == "P" & Genotype == "110" ~ 4,
    Species == "P" & Genotype == "101" ~ 5,
    Species == "P" & Genotype == "011" ~ 6,
    Species == "P" & Genotype == "111" ~ 7,
    TRUE ~ -999
  ))

tmin <- 90000
tmax <- 92000

p_full <- ggplot(daten_longer_cumsum, 
                 aes(x=time)) + 
  geom_ribbon(aes(ymin = lower_freq, ymax = cumsum_freq, fill= Genotype), alpha = 1) + 
  facet_grid(rows =vars(Species_label)) + 
  scale_fill_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  theme_minimal() +
  ylab("cumulative frequency") +
  ggtitle("time span: 0 - 1,000") +
  theme(strip.text.y = element_text(size = 10, face = 2, hjust =0),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size= 12),
        legend.title = element_text(size = 13, face =2)
  ) + 
  guides(fill= guide_legend(override.aes = list(size = 10))) 


p1 <- ggplot(daten_longer_cumsum |> 
               filter(time <=1000), 
       aes(x=time)) + 
  geom_ribbon(aes(ymin = lower_freq, ymax = cumsum_freq, fill= Genotype), alpha = 1 ) + 
  facet_grid(rows =vars(Species_label)) + 
  scale_fill_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  theme_minimal() +
  ylab("cumulative frequency") +
  ggtitle("time span: 0 - 1,000") +
  theme(strip.text.y = element_text(size = 10, face = 2, hjust =0),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size= 12),
        legend.title = element_text(size = 13, face =2)
        ) + 
  guides(fill= guide_legend(override.aes = list(size = 10))) 

p2plt <- TRUE

if(any(daten_longer_cumsum[["time"]] > 9000 & daten_longer_cumsum[["time"]] < 10000 ) ){
  p2 <- ggplot(daten_longer_cumsum |> filter(time > 9000 & time < 10000), 
               aes(x=time)) + 
    geom_ribbon(aes(ymin = lower_freq, ymax = cumsum_freq, fill= Genotype), alpha = 1 ) + 
    facet_grid(rows =vars(Species_label)) + 
    scale_fill_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
    theme_minimal() +
    ylab("cum. freq") +
    ggtitle("time span: 9,000 - 10,000") +
    ylab("cumulative frequency") +
    theme(strip.text.y = element_text(size = 10, face = 2, hjust = 0),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.text = element_text(size= 12),
          legend.title = element_text(size = 13, face =2)
    ) +
    guides(fill= guide_legend(override.aes = list(size = 10)))#+
}else{
  p2 <- ggplot(daten_longer_cumsum |> filter(time > 1000), 
               aes(x=time)) + 
    geom_ribbon(aes(ymin = lower_freq, ymax = cumsum_freq, fill= Genotype), alpha = 1 ) + 
    facet_grid(rows =vars(Species_label)) + 
    scale_fill_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
    theme_minimal() +
    ylab("cum. freq") +
    ggtitle("time span: 1000 - 100,000") +
    ylab("cumulative frequency")
    theme(strip.text.y = element_text(size = 10, face = 2, hjust = 0),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.text = element_text(size= 12),
          legend.title = element_text(size = 13, face =2)
    )
    guides(fill= guide_legend(override.aes = list(size = 10))) #+
}

p3plt <- FALSE

if(any(daten_longer_cumsum[["time"]] > 99000 & daten_longer_cumsum[["time"]] < 100000 ) ){
  p3plt <- TRUE
  p3 <- ggplot(daten_longer_cumsum |> filter(time > 99000 & time < 100000), 
               aes(x=time)) + 
    geom_ribbon(aes(ymin = lower_freq, ymax = cumsum_freq, fill= Genotype), alpha = 1 ) + 
    facet_grid(rows =vars(Species_label)) + 
    scale_fill_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
    theme_minimal() +
    ylab("cum. freq") +
    ggtitle("time span: 99,000 - 100,000") +
    ylab("cumulative frequency") +
    theme(strip.text.y = element_text(size = 10, face = 2, hjust = 0),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.text = element_text(size= 12),
          legend.title = element_text(size = 13, face =2)#,
    ) + 
    guides(fill= guide_legend(override.aes = list(size = 10)))#+
}




if(!p2plt & !p3plt){out <- p1 + plot_layout(nrow = 3, guides = 'collect')}
if(!p2plt & p3plt){out <- p1 + p3  + plot_layout(nrow = 3, guides = 'collect')}
if(p2plt & !p3plt){out <- p1 + p2  + plot_layout(nrow = 3, guides = 'collect')}
if(p2plt & p3plt){out <- p1 + p2 + p3 + plot_layout(nrow = 3, guides = 'collect')}
out

svg("testout.svg", width = 8, height = 12)
print(out)
dev.off()

if(exists("p1") & exists("p2") & exists("p3") & exists("dynamicflt_plt1") & exists("dynamicflt_plt1") & exists("dynamicflt_plt3")){
  Fig1_like <- (dynamicflt_plt1 + p1) / (dynamicflt_plt2 + p2)+ (dynamicflt_plt3 + p3) + 
    plot_layout(guides = "collect") #+ plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
  
  svg('outextended.svg', width = 10, height = 10)
  print(Fig1_like)
  dev.off()
}

# Generate Figure 1 
if(CH1 == 0.1 & CP1 == 0.1 & CP2 == 3 & CH2 == 3 & phi == 0.5 & seed == 1600){
  Fig1 <- (dynamicflt_plt1 + p1) / (dynamicflt_plt2 + p2)+ (dynamicflt_plt3 + p3) + 
    plot_layout(guides = "collect") + plot_annotation(tag_levels = "i", tag_prefix = "(", tag_suffix = ")")
  
  pdf(paste0(plotdir,"/Fig1.pdf"), width = 10, height = 10)
  print(Fig1)
  dev.off()
  
}
