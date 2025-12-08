suppressPackageStartupMessages(require("tidyverse"))
suppressPackageStartupMessages(require("paletteer"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("patchwork"))
suppressPackageStartupMessages(require("ggpubr"))



dynamics_fct <- function(daten_longer_cumsum, tmin, tmax){
  plt <- ggplot(daten_longer_cumsum |>
         filter(time < tmax & time > tmin & freq > 0),
       aes(x=time)) +
  geom_line(aes(x=time, y=freq, linetype = Species_label, color = Genotype, alpha = Species_label)) +
  scale_alpha_manual(values=c(1, 0.5, 1)) +
  geom_point(aes(x=time, y=freq, shape = Species_label, color = Genotype, alpha = Species_label), size = 0.4) +
  facet_grid(rows = vars(IPart), scales = "free_y") +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  theme_bw() +
  ylab("freq") +
  scale_linetype_manual(values = c("Host H"="solid", "Host M"="11", "Pathogen"="solid")) +
  ylab("Frequency") +
  ggtitle(paste0("time span: ",format(tmin, big.mark = ",", scientific = FALSE)," - ",format(tmax, big.mark = ",", scientific = FALSE))) +
  ylim(0,1)  +
  theme() +
  scale_x_continuous(limits = c(tmin, tmax), expan = expansion())   +
  theme(plot.margin = unit(c(0.05, 0.2, 0, 1), 'lines'), 
        panel.grid.major.x = element_line(color = "gray40"), 
        panel.grid.minor.x = element_line(color = "gray80"),
        legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(hjust = 0, size = 12, face = 2),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))+ 
    guides(color = "none", alpha = "none", linetype = guide_legend(override.aes = list(color = "gray50", size = c(4,4,4), linewidth = c(2,1,2), alpha = c(1,0.8,1))))
  return(plt)
}

cumdynamics_plt <- function(daten_longer_cumsum, tmin, tmax){
  plt <- ggplot(daten_longer_cumsum |> 
           filter(time <= tmax & time >= tmin), 
         aes(x=time)) + 
    geom_ribbon(aes(ymin = lower_freq, ymax = cumsum_freq, fill= Genotype), alpha = 1 ) + 
    facet_grid(rows =vars(Species_label)) + 
    scale_fill_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                  "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
    theme_minimal() +
    ylab("cumulative frequency") +
    ggtitle(paste0("time span: ",format(tmin, big.mark = ",", scientific = FALSE)," - ",format(tmax, big.mark = ",", scientific = FALSE))) +
    theme(strip.text.y = element_text(size = 10, face = 2, hjust =0),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.text = element_text(size= 12),
          legend.title = element_text(size = 13, face =2)
    ) + 
    guides(fill= guide_legend(override.aes = list(size = 10))) 
  return(plt)
}



# Script to generate Fig 1D and 1E

prefix <- "ext_reps"
suffix <- "2025-01-07"
omegaH <- 0.1
omegaP <- 0.1
XiH <- 3
XiP <- 3
phi <- 0.5
seed <- 1600
plotdir <- "Figures/main_figures"


outplt <- paste0(plotdir,"/dynamics_cH2_",XiH,"_cP2_",XiP,"_phi_",phi,"_cH1_",omegaH,"_cP1_",omegaP,".pdf")
outplt_extended <- paste0(plotdir,"/dynamics_cH2_",XiH,"_cP2_",XiP,"_phi_",phi,"_cH1_",omegaH,"_cP1_",omegaP,"_extended.pdf")

if(!dir.exists(plotdir)){
  dir.create(plotdir, recursive = TRUE)
}

# These parameters are fixed throughout the simulations
sigma <- 0.85
betaH <- 1
betaP <- 1


# Simulation prefix. This is the entire file name prefix for all the result files of the simulation
simprefix <- paste0(prefix,"_phi_",phi,"_cH_1_",omegaH,"_cH_2_",XiH,"_cP_1_",omegaP,"_cP_2_",XiP,"_s_",sigma,"_bH_",betaH,"_bP_",betaP,"_seed_",seed,"_",suffix)

# Read the dynamics from the file
daten_longer_cumsum <- read_tsv(paste0("Data_Fig1/",simprefix,"_long.tsv"))

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


#######################
# Fig 1E ##############
#######################

dynamicflt_plt_full <- ggplot(daten_longer_cumsum, aes(x=time)) +
  # add grey rectangles to indicate the time intervals for which the dynamics are being plotted in Fig 1F
  annotate("rect", fill = "gray20", alpha = 0.5, xmin = c(0, 9000, 99000), xmax = c(1000, 10000, 100000), ymin = -Inf, ymax = Inf) +  
  # Plot the genotype dynamics as line
  geom_line(aes(x=time, y=freq, color = Genotype),  linewidth = 0.2) +
  # Generate a facet grid with species on rows
  facet_grid(rows = vars(Species_label)) +
  # Define the color scale
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  # Set the base size of all plot text to 14
  theme_bw(base_size = 14) +
  # Set the axis labels
  ylab("freq") +
  ylab("Genotype frequency") +
  ylim(0,1)  +
  scale_x_continuous(limits = c(0, 100000), expan = expansion())   +
  # Fine tune the appearance of the plot
  theme(plot.margin = unit(c(0.05, 0.2, 0, 1), 'lines'),
        panel.grid.minor.x = element_line(color = "gray80"),
        legend.key.width = unit(1,"cm"),
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
  
svg(paste0(plotdir,"/fig1_fulldynamics.svg"), width = 14, height = 4)
print(dynamicflt_plt_full)
dev.off()


#####################################################
# Fig 1F fluctuations of allele frequencies 0-1,000 #
#####################################################

dynamicflt_plt1 <- dynamics_fct(daten_longer_cumsum, 0, 1000)
dynamicflt_plt2 <- dynamics_fct(daten_longer_cumsum, 9000, 10000)
dynamicflt_plt3 <- dynamics_fct(daten_longer_cumsum, 99000, 100000)


# dynamicflt_plt1 <- ggplot(daten_longer_cumsum |>
#                            filter(time < 1000 & freq > 0),
#                          aes(x=time)) +
#   geom_line(aes(x=time, y=freq, linetype = Species_label, color = Genotype, alpha = Species_label)) +
#   scale_alpha_manual(values=c(1, 0.5, 1)) +
#   geom_point(aes(x=time, y=freq, shape = Species_label, color = Genotype, alpha = Species_label), size = 0.4) +
#   facet_grid(rows = vars(IPart), scales = "free_y") +
#   scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
#                                  "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
#   theme_bw() +
#   ylab("freq") +
#   scale_linetype_manual(values = c("Host H"="solid", "Host M"="11", "Pathogen"="solid")) +
#   ylab("Frequency") +
#   ggtitle("time span: 0 - 1,000") +
#   ylim(0,1)  +
#   theme() +
#   scale_x_continuous(limits = c(0, 1000), expan = expansion())   +
#   theme(plot.margin = unit(c(0.05, 0.2, 0, 1), 'lines'), 
#         panel.grid.major.x = element_line(color = "gray40"), 
#         panel.grid.minor.x = element_line(color = "gray80"),
#         legend.key.width = unit(1,"cm"),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         strip.background = element_blank(), 
#         strip.text = element_text(hjust = 0, size = 12, face = 2),
#         axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
#         axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 12))+ 
#   guides(color = "none", linetype = "none", alpha = "none", shape = "none")
# 
# 
# dynamicflt_plt2 <- ggplot(daten_longer_cumsum |>
#          filter(time < 10000 & time > 9000 & freq > 0),
#        aes(x=time)) +
#   geom_line(aes(x=time, y=freq, linetype = Species_label, color = Genotype, alpha = Species_label)) +
#   scale_alpha_manual(values=c(1, 0.5, 1)) +
#   #geom_line(aes(x=time, y=freq, color = Genotype, group = interaction(Species, Genotype))) +
#   geom_point(aes(x=time, y=freq, shape = Species_label, color = Genotype, alpha = Species_label), size = 0.4) +
#   #geom_line(aes(y = cumsum_freq, color = Genotype)) +
#   facet_grid(rows = vars(IPart), scales = "free_y") +
#   scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
#                                  "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
#   theme_bw() +
#   theme(panel.grid.major.x = element_line(color = "gray20"), panel.grid.minor.x = element_line(color = "gray80")) +
#   ylab("freq") +
#   #labs(subtitle = "(b)") +
#   scale_linetype_manual(values = c("Host H"="solid", "Host M"="11", "Pathogen"="solid")) +
#   ylab("Frequency") +
#   ggtitle("time span: 9,000 - 10,000") +
#   ylim(0,1) +
#   theme(plot.margin = unit(c(0.05,0.2,0,1), 'lines')) + 
#   scale_x_continuous(limits = c(9000, 10000), expan = expansion())    +
#   theme(legend.key.width = unit(1,"cm"),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         strip.background = element_blank(), 
#         strip.text = element_text(hjust = 0, size = 12, face = 2),
#         axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
#         axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 12)) +
#   guides(color = "none", alpha = "none", 
#          linetype = guide_legend(override.aes = list(color = "gray50", size = c(4,4,4), linewidth = c(2,1,2), alpha = c(1,0.8,1))))
# 
# 
# dynamicflt_plt3 <- ggplot(daten_longer_cumsum |>
#                            filter(time > 99000 & time < 100000 & freq > 0),
#                          aes(x=time)) +
#   geom_line(aes(x=time, y=freq, linetype = Species_label, color = Genotype, alpha = Species_label)) +
#   scale_alpha_manual(values=c(1, 0.5, 1)) +
#   #geom_line(aes(x=time, y=freq, color = Genotype, group = interaction(Species, Genotype))) +
#   geom_point(aes(x=time, y=freq, shape = Species_label, color = Genotype, alpha = Species_label), size = 0.4) +
#   #geom_line(aes(y = cumsum_freq, color = Genotype)) +
#   facet_grid(rows = vars(IPart), scales = "free_y") +
#   scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
#                                  "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
#   theme_bw() +
#   theme(panel.grid.major.x = element_line(color = "gray40"), panel.grid.minor.x = element_line(color = "gray80")) +
#   ylab("freq") + #+
#   #labs(subtitle = "(c)") +
#   scale_linetype_manual(values = c("Host H"="solid", "Host M"="11", "Pathogen"="solid")) +
#   ylab("Frequency") +
#   ylim(0,1) +
#   ggtitle("time span: 99,000 - 100,000") +
#   theme(plot.margin = unit(c(0.05, 0.2,0,1), 'lines')) +
#   scale_x_continuous(limits = c(99000, 100000), expan = expansion()) +
#   theme(legend.key.width = unit(1,"cm"),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(),
#         strip.background = element_blank(), 
#         strip.text = element_text(hjust = 0, size = 12, face = 2),
#         axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
#         axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 12)) + 
#   guides(color = "none", linetype = "none", alpha = "none", shape = "none")


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


p1 <- cumdynamics_plt(daten_longer_cumsum, tmin = 0 , tmax = 1000)
p2 <-  cumdynamics_plt(daten_longer_cumsum, tmin = 9000, tmax = 10000)
p3 <-  cumdynamics_plt(daten_longer_cumsum, tmin = 99000, tmax = 100000)

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





Fig1 <- (dynamicflt_plt1 + p1) / (dynamicflt_plt2 + p2)+ (dynamicflt_plt3 + p3) + 
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "i", tag_prefix = "(", tag_suffix = ")")
  
svg(paste0(plotdir,"/Fig1.svg"), width = 10, height = 10)
print(Fig1)
dev.off()
  

