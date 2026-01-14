#require("grid") # for doing some grid work
#require("gridExtra") # for simple ggplot layout tasks
require("patchwork") # needed for better plot layouts
require("egg") # needed for tag_facet
require("ggExtra") # for the ggMarginal function
require("ggpubr") # For annotations

# Description: 
# This script will plot the proportion of simulations for which genotype with given numbers of resistance and or virulence
# alleles are maintained
# For example: 0 + 1 indicates that the haplotype wi

# Directory with the data files for plotting
indir <- "results_random"

# Results for the hosts and the pathogen
datH_dat <- read_tsv(paste0(indir, "/datH_dat.tsv"))
datP_dat <- read_tsv(paste0(indir, "/datP_dat.tsv"))

# Set output directory (note)
outdir <- "Figures"

main_figures <- paste0(outdir,"/main_figures")
supplementary_figures <- paste0(outdir, "/supplementary_figures")
tabledir <- "supplementary_tables"

# Check if all output directories already exist. If not, create them
# parent output directory
if(!dir.exists(outdir)){
  dir.create(outdir)
}

# directory for the main figures
if(!dir.exists(main_figures)){
  dir.create(main_figures)
}

# directory for the supplementary figures
if(!dir.exists(supplementary_figures)){
  dir.create(supplementary_figures)
}

# directory for the supplementary tables
if(!dir.exists(tabledir)){
  dir.create(tabledir)
}

datH_dat <- datH_dat |>
  # Convert to factor to ensure consistent plotting
  mutate(combi_maintained = factor(combi_maintained, levels = 
                                     c("private: S, shared: S",
                                       "private: R, shared: S",
                                       "private: R/S, shared: S",
                                       "private: S, shared: R/S",                                       
                                       "private: S, shared: R",
                                       "private: R, shared: R",
                                       "private: R/S, shared: R",                                      
                                       "private: R, shared: R/S",
                                       "private: R/S, shared: R/S")))


###################################
# Generate subfigures S4 and S4b ##
###################################

XiH_value <- 3
XiP_value <- 3

caseH_phi <- datH_dat |>
  filter(XiH == XiH_value & XiP == XiP_value) |>
  mutate(phi = factor(phi, levels=c(0.5,0.7,0.9,1)))


caseP_phi <- datP_dat |>
  filter(XiH == XiH_value & XiP == XiP_value) |>
  mutate(phi = factor(phi, levels=c(0.5,0.7,0.9,1)))


# Summarize the data for annotating the barplot in (b)
maintained_summary <- caseH_phi |> 
  group_by(phi, Species_label, combi_maintained) |>
  tally()  |>
  group_by(phi, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))



############################
## Plot S4Fig A ###########
###########################

maintained_alleles_standard <- ggplot(caseH_phi, aes(x=OmegaH, y=OmegaP)) + 
  geom_point(aes(color = combi_maintained), show.legend=T, size = 0.05) + 
  facet_grid(cols = vars(Species_label),
             rows = vars(phi),
             labeller=label_bquote(rows = phi == .(as.character(phi)))) +
  theme_bw() + 
  scale_color_manual(values = c("gray40", "mediumvioletred","#925E9FFF","mediumseagreen","#FDAF91FF","#AD002AFF","cadetblue1","#0099B4FF","#00468BFF","lightyellow"),
                      name = "maintained alleles\nin population",drop =FALSE) +
  guides(color=guide_legend(override.aes = list(size = 5, shape = 15))) +
  labs(x = bquote(bold("maximum cost of resistance "*Omega[H])), 
       y = bquote(bold("maximum cost of virulence "*Omega[P]))) +
  # Modify the appearnce of the plot
  theme_bw() +
  theme(strip.text = element_text(size = 16, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = 2),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(vjust = -0.75, size = 15),
        axis.title.y = element_text(size = 15),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.spacing = unit(0.5, "cm"),
        legend.position = "none")  +
  # Set the axis 
  scale_x_continuous(limits = c(0, 0.3),
                     labels = c(0, 0.1, 0.2, 0.3),
                     expand = expansion(c(0,0), c(0,0.005)),
                     breaks = c(0, 0.1, 0.2, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3),
                     labels = c(0, 0.1, 0.2, 0.3),
                     expand = expansion(c(0,0), c(0,0.005)),
                     breaks = c(0, 0.1, 0.2, 0.3)) 


#############################
# Generate S4Fig B #########
############################

f_combi <- ggplot(caseH_phi , aes(x="", fill=combi_maintained)) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion of 10,000 simulations") +
  theme_bw() +
  facet_grid(rows = vars(phi), cols = vars(Species_label), labeller=label_bquote(rows = phi == .(as.character(phi)) , cols = .(Species_label))) + 
  scale_fill_manual( values = c("gray40", "mediumvioletred","#925E9FFF","mediumseagreen","#FDAF91FF","#AD002AFF","cadetblue1","#0099B4FF","#00468BFF","lightyellow"),name = "maintained alleles\nin population",drop =FALSE) +
  labs(title = "", 
       x = "", 
       y = "Proportion out of\n10,000 simulations")+
  guides(fill="none") + 
  geom_text(data = maintained_summary ,aes(y=prop_y, label = label, group = combi_maintained), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 5,
            color = "white") +
  theme(plot.margin = margin(0,0,0,0, 'cm'),
        strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")), 
        axis.text = element_text(size = 14),
        axis.title.x = element_text(vjust = -0.75, size = 15),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm"), size = 15),
        panel.grid = element_blank(),
        panel.border = element_blank(),  # Removes the box around the facets
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        title = element_text(face = "bold", size = 16),
        panel.spacing.x = unit(0, "cm"),
        panel.spacing.y = unit(0.5, "cm"),
        panel.background = element_rect(fill = "gray60"),
        legend.text = element_text(size = 14))




# Summary for the maintained alleles
combi_summary <- datH_dat %>% 
  group_by(phi_label, Species_label, combi_maintained) %>% 
  tally()  %>%
  group_by(phi_label, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))


#############################
# Generate S4Fig C #########
############################

combimaintained_acrossall_plt <- ggplot(datH_dat, aes(x="", fill= factor(combi_maintained))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(rows = vars(phi_label), cols =vars(Species_label), labeller=labeller(phi_label =label_parsed)) +
  scale_fill_manual(
    values = c("gray40", "mediumvioletred","#925E9FFF","mediumseagreen","#FDAF91FF","#AD002AFF","cadetblue1","#0099B4FF","#00468BFF","lightyellow"),
    name = "maintained alleles\nin population",
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(x = "", 
       y = "Proportion of\n90,000 simulations") +
  theme_bw() +
  geom_text(data = combi_summary, aes(y=prop_y, label = label, group = combi_maintained), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 5,
            color = "white") +
  theme(strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")), 
        axis.text = element_text(size = 14),
        axis.title.x = element_text(vjust = -0.75, size = 15),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm"), size = 15),
        panel.grid = element_blank(),
        panel.border = element_blank(),  # Removes the box around the facets
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        title = element_text(face = "bold", size = 16),
        panel.spacing.y = unit(0.5, "cm"),
        panel.spacing.x = unit(0, "cm"),
        panel.background = element_rect(fill = "gray60")) +
  #strip.background = element_rect(fill = "white", colour = "white")) +
  ggtitle("Across cost function shapes\nand maximum costs") +
  theme(legend.text = element_text(size = 13)) 

## Combine all plots into the final figure
SFig4_new <- maintained_alleles_standard + theme(legend.position = "none") + ggtitle("Effect of maximum costs\nfor fixed shapes") +
  plot_spacer() +
  f_combi + theme(legend.position = "none") +  ggtitle("Summary maximum\ncosts fixed shapes")  +
  combimaintained_acrossall_plt +
  plot_layout(
    widths = c(2, 0.01,  1, 1),    # 1/3 : 2/3 column ratio
    heights = c(1)    # equal row heights
  ) +
  plot_annotation(
    tag_levels = list(c("(a)", "(b)", "(c)")),
    tag_prefix = "",
    tag_suffix = ""
  ) &
  theme(
    plot.tag = element_text(face = 2, size = 14),
    plot.title = element_text(face = 2, size = 16)
  )

pdf(paste0(supplementary_figures,"/S4Fig_alleles_maintained.pdf"), width = 15, height = 9)
print(SFig4_new)
dev.off()

ggsave(filename = paste0(supplementary_figures,"/S4Fig_alleles_maintained.png"), 
       plot = SFig4_new,
       width = 15, 
       height = 9, 
       units = "in", dpi = 400)

#################################
## Figure S5
#################################


# Now make a more detailed version of the plot. First generate the corresponding summary
# Note set the label threshold to 5% here
combi_summary_precise <- datH_dat |>
  mutate(phi = factor(phi)) |>
  group_by(XiH, XiP, phi, Species_label, combi_maintained) |> 
  tally()  |>
  ungroup() |>
  group_by(XiH, XiP, phi, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))


################################
## Plot S5Fig A ###############
###############################

Rcombi_costcombi_plt_hostH <- ggplot(datH_dat |> filter(Species_label == "Host H"), aes(x=as.factor(phi), 
                                                                                        fill= factor(combi_maintained))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Proportion of host H", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols = vars(XiH),rows = vars(XiP), labeller=label_bquote(cols = Xi[H] == .(XiH), rows = Xi[P] == .(XiP))) +
  scale_fill_manual(
    values = c("gray40", "mediumvioletred","#925E9FFF","mediumseagreen","#FDAF91FF","#AD002AFF","cadetblue1","#0099B4FF","#00468BFF","lightyellow"),
    name = "maintained alleles\nin population",
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(x = bquote("Relative proportion of host species H ("*phi*")"), 
       y = "Proportion out of\n10,000 simulations") +
  theme_bw() +
  geom_text(data = combi_summary_precise |> filter(Species_label == "Host H"), aes(y=prop_y, label = label, group = combi_maintained), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 4,
            color = "white") +
  theme(strip.text = element_text(size = 16, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        legend.text = element_text(size = 12),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_text(size = 14, face = 2),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -0.75, size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        title = element_text(size = 16, face = 2)) +
  ggtitle("Host species H")

################################
## Plot S5Fig B ###############
###############################

Rcombi_costcombi_plt_hostM <- ggplot(datH_dat |> mutate(phi = factor(phi)) |> filter(Species_label == "Host M"),
                                     aes(x=phi, fill= factor(combi_maintained))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Proportion of host H", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols = vars(XiH),rows = vars(XiP), labeller=label_bquote(cols = Xi[H] == .(XiH), rows = Xi[P] == .(XiP))) +
  scale_fill_manual(
    values = c("gray40", "mediumvioletred","#925E9FFF","mediumseagreen","#FDAF91FF","#AD002AFF","cadetblue1","#0099B4FF","#00468BFF","lightyellow"),
    name = "maintained alleles\nin population",
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(x = bquote("Relative proportion of host species H ("*phi*")"), 
       y = "Proportion out of\n10,000 simulations") +
  theme_bw() +
  geom_text(data = combi_summary_precise |> filter(Species_label == "Host M"), aes(y=prop_y, label = label, group = combi_maintained), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 4,
            color = "white") +
  theme(strip.text = element_text(size = 16, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        legend.text = element_text(size = 12),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_text(size = 14, face = 2),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -0.75, size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        title = element_text(size = 16, face = 2)) +
  ggtitle("Host species M")#+

#################################
## Generate the combined plot ###
#################################

SFig5_new <- Rcombi_costcombi_plt_hostH + theme(legend.position = "none") +
  Rcombi_costcombi_plt_hostM  +
  plot_layout(
    widths = c(1, 1),    # 1/3 : 2/3 column ratio
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = list(c("(a)", "(b)")),
    tag_prefix = "",
    tag_suffix = ""
  ) &
  theme(
    plot.tag = element_text(face = 2, size = 14),
    plot.title = element_text(face = 2, size = 16)
  )

pdf(paste0(supplementary_figures,"/S5Fig_maintained_detail.pdf"),width = 15, height=8)
print(SFig5_new)
dev.off()

ggsave(filename = paste0(supplementary_figures,"/S5Fig_maintained_detail.png"), 
       plot = SFig5_new,
       width = 15, 
       height = 8, 
       units = "in", dpi = 400)

