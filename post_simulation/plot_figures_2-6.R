require("tidyverse")
require("argparser")
require("ggpubr") # for using ggtext, gradient_fill, gradient_color etc. and to add lm results etc
require("grid") # for doing some grid work
require("paletteer") # for more colors
require("gridExtra") # for simple ggplot layout tasks
require("ggridges") # if ridges to be plotted
require("viridis") # needed for access to color blind friendly palettes
require("hrbrthemes") # needed for access to additional themes
require("patchwork") # needed for better plot layouts
require("egg") # needed for tag_facet
require("ggExtra")


args <- arg_parser("Command line argument parser")
args <- add_argument(args, arg="--parentdir",
                     help="Suffix to append to all filenames",
                     type="character",
                     default = getwd())
args <- add_argument(args, arg="--suffix",
                     help="Suffix to append to all filenames",
                     type="character",
                     default= "2025-01-27")
args <- add_argument(args, arg="--outdir",
                     help="output directory for the plots. Relative to parent dir.",
                     default = "../Figures",
                     type="character")
args <- add_argument(args, arg="--tabledir",
                     help="output directory for the supplementary tables. Relative to parent dir.",
                     default = "../supplementary_tables",
                     type="character")
pargs <- parse_args(args)

# Parse the command line arguments
indir <- pargs[["parentdir"]]
outdir <- pargs[["outdir"]]
suffix <- pargs[["suffix"]]
tabledir <- pargs[["tabledir"]]


setwd(indir)
getwd()

main_figures <- paste0(outdir,"/main_figures")
supplementary_figures <- paste0(outdir, "/supplementary_figures")

# Check if all output directories already exist. If no, create them
if(!dir.exists(outdir)){
  dir.create(outdir)
}

if(!dir.exists(main_figures)){
  dir.create(main_figures)
}

if(!dir.exists(supplementary_figures)){
  dir.create(supplementary_figures)
}

if(!dir.exists(tabledir)){
  dir.create(tabledir)
}

# Create an overview of all files in end summaries
files_toread <- list.files("end_summaries", full.names =T)

# Check which files to read for the host side
to_readH <- files_toread[grep(glob2rx("outH*.tsv"), list.files("end_summaries"))]
# Check which files to read for the pathogen side
to_readP <- files_toread[grep(glob2rx("outP*.tsv"), list.files("end_summaries"))]

# Read all host results
datH <- lapply(to_readH, function(x){
  out <- read_tsv(x)
})

# Combine all host results
datH_dat <- bind_rows(datH)

head(datH_dat)
names(datH_dat)

# Add a couple of new columns to the host master file
datH_dat <- datH_dat |>
  mutate(shape_combi=paste(CH2, CP2, sep="_")) |>
  mutate(combi_label = paste("atop(c[H]^{(2)} ==", CH2,",", "c[P]^{(2)} == ", CP2,")")) |>
  mutate(phi_label = paste("phi ==", phi)) |>
  mutate(Species_label = case_when(
    Species == "H1" ~ "Host H",
    Species == "H2" ~ "Host M",
    TRUE ~ "other"
  ))

datH_dat <- datH_dat |>
  mutate(combi_maintained = case_when(
    R_alleles == 0 ~ "private: S, shared: S",
    R_alleles == 1 & poly == 0 ~ "private: R, shared: S",
    R_alleles == 1 & poly == 1 ~ "private: R/S, shared: S",
    R_alleles == 2 & poly == 0 ~ "private: S, shared: R",
    R_alleles == 2 & poly == 2 ~ "private: S, shared: R/S",
    R_alleles == 3 & poly == 0 ~ "private: R, shared: R",
    R_alleles == 3 & poly == 1 ~ "private: R/S, shared: R",
    R_alleles == 3 & poly == 2 ~ "private: R, shared: R/S",
    R_alleles == 3 & poly == 3 ~ "private: R/S, shared: R/S",
    TRUE ~ "forgotten"
  )) |>
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

# Read the file for the pathogen
datP <- lapply(to_readP, function(x){
  out <- read_tsv(x)
})

# Combine all files for the pathogen
datP_dat <- bind_rows(datP)
# Add more information to the pathogen df
datP_dat <- datP_dat |>
  mutate(shape_combi=paste(CH2, CP2, sep="_")) |>
  mutate(Species_label="Pathogen")

# Check if everything looks right
# There should be a total of 630,000 lines in the combined host data frame
# Explanation: 
# 9 (cH2 and cP2) combinations x 10,000 (random cH1, cP1 combinations) x 3 (phi values: 0.5, 0.7, 0.9) x 2 hosts = 540,000
# 9 (cH2 and cP2) combinations x 10,000 (random cH1, cP1 combinations) x 1 (phi values:1 ) x 1 host = 90,000
# Total: 630,000
# There should be a total of 360,000 lines in the pathogen combined data frame
# 9 (cH2 and cP2) combinations x 10,000 (random cH1, cP1 combinations) x 4 (phi values: 0.5, 0.7, 0.9, 1.0) = 360,000
nrow(datH_dat)
nrow(datP_dat)
datH_dat |> 
  group_by(shape_combi, phi, Species_label) |>
  summarize(n=n()) |> filter(n !=10000)

datP_dat |> 
  group_by(shape_combi, phi, Species_label) |>
  summarize(n=n()) |> filter(n !=10000)



# Figure 2 ----------------------------------------------------------------

## Analyze the standard case

standard_caseH <- datH_dat |>
  filter(CH2==3 & CP2==3 & phi == 0.5) 
standard_caseH <- standard_caseH |>
  mutate(maxR = sapply(standard_caseH[["range_summary"]], function(x){max(which(intToBits(x)==1))-1}))

standard_caseP <- datP_dat |>
  filter(CH2==3 & CP2==3 & phi == 0.5) 

standard_caseP <- standard_caseP |>
  mutate(maxV = sapply(standard_caseP[["range_summary"]], function(x){max(which(intToBits(x)==1))-1}))



rdens_std_plt <- ggplot(standard_caseH,  aes(x = CH1, color = as.factor(maxR), fill = as.factor(maxR))) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values= c("0" = "#9EBCDA",
                              "1" = "#8C6BB1",
                              "2" = "#810F7C",
                              "3" = "#4D004B"), name = "maximum number of\nR alleles/haplotype") +
  scale_color_manual(values= c("0" = "#9EBCDA",
                               "1" = "#8C6BB1",
                               "2" = "#810F7C",
                               "3" = "#4D004B"), name = "maximum number of\nR alleles/haplotype")  +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 14)) + 
  ylab("Probability\ndensity") +
  theme(plot.margin = margin(0,0,0,0, 'cm')) +
  xlab(bquote("Maximum cost of resistance"~Omega[H]))


rdensCP_std_plt <- ggplot(standard_caseH,  aes(x = CP1, color = as.factor(maxR), fill = as.factor(maxR))) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values= c("0" = "#9EBCDA",
                              "1" = "#8C6BB1",
                              "2" = "#810F7C",
                              "3" = "#4D004B"), name = "maximum number of\nR alleles/haplotype") +
  scale_color_manual(values= c("0" = "#9EBCDA",
                               "1" = "#8C6BB1",
                               "2" = "#810F7C",
                               "3" = "#4D004B"), name = "maximum number of\nR alleles/haplotype")  +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_text(size = 16), axis.text = element_text(size = 14)) +
  theme(legend.position = "none") + 
  ylab("Probability\ndensity") +
  theme(plot.margin = margin(0,0,0,0, 'cm')) +
  xlab(bquote("Maximum cost of virulence"~Omega[P])) #+



rpnt_std_plt <- ggplot(standard_caseH,  aes(x = CH1, CP1, color = as.factor(maxR))) +
  geom_point(alpha=0.6, size = 0.8) +
  scale_color_manual(values= c("0" = "#9EBCDA",
                               "1" = "#8C6BB1",
                               "2" = "#810F7C",
                               "3" = "#4D004B"), name = "maximum number of\nR alleles/haplotype")  +
  theme_bw(base_size = 14) +
  ylab(bquote("Maximum cost of virulence"~Omega[P])) +
  xlab(bquote("Maximum cost of resistance"~Omega[H])) +
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 15))) +
  theme( axis.text = element_text(size = 14),
         axis.title.x = element_text(vjust = -0.75, size = 16),
         axis.title.y = element_text(size = 16))#+


outgg <- ggMarginal(rpnt_std_plt + 
                      scale_color_manual(values= c("0" = "#9EBCDA",
                                                   "1" = "#8C6BB1",
                                                   "2" = "#810F7C",
                                                   "3" = "#4D004B"), name = "maximum number of R alleles/haplotype") +
                      theme(axis.title = element_text(size = 16), legend.title = element_text(size = 18),
                            legend.text = element_text(size = 18)), type="density", groupColour = TRUE, groupFill = TRUE) 


pdf(paste0(main_figures,"/Fig_2b_std_case_rpnt_CH2_3_CP2_3_phi_0.5.pdf"), width = 7, height = 6.5)
print(outgg)
dev.off()

composite_nb_alleles  <- (rdens_std_plt + scale_x_continuous(position = "top")) + plot_spacer() + rpnt_std_plt + (rdensCP_std_plt + 
                                                                                                                    scale_y_continuous(position = "right") + 
                                                                                                                    scale_x_continuous(position = "top") + rotate()) +
  plot_layout(widths = c(4, 1), heights = c(1,4))  & 
  plot_annotation(tag_levels = list(c("(a)","(b)","(c)",""))) &
  theme(plot.tag = element_text(face = 2))

pdf(paste0(main_figures,"/Supp_Fig_3_std_case_composite_CH2_3_CP2_3_phi_0.5.pdf"), width = 8.5, height = 8.5)
print(composite_nb_alleles )
dev.off()



# Figure 2 ----------------------------------------------------------------

CH2_value <- 3
CP2_value <- 3
    
caseH_phi <- datH_dat |>
    filter(CH2 == CH2_value & CP2== CP2_value) |>
    mutate(phi = factor(phi, levels=c(0.5,0.7,0.9,1)))
    
    
caseP_phi <- datP_dat |>
  filter(CH2== CH2_value & CP2 == CP2_value) |>
  mutate(phi = factor(phi, levels=c(0.5,0.7,0.9,1)))
    
    
Ralleles_summary <- caseH_phi %>% 
  group_by(phi, Species_label, R_alleles = factor(R_alleles, levels = rev(as.character(0:3)))) %>% 
  tally()  %>%
  group_by(phi, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
   mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  )) |>
  mutate(Rcombi = case_when(
    R_alleles == 0 ~ "none",
    R_alleles == 1 ~ "private only",
    R_alleles == 2 ~ "ancestral only",
    R_alleles == 3 ~ "both", 
    TRUE ~ "unaccounted"
  ))

out_std_Rsum <- Ralleles_summary |>
  select(Species_label, phi, Rcombi, prop_y) |>
  mutate(prop_y = format(round(prop_y, digits =3), nsmall = 3)) |>
  pivot_wider(values_from = prop_y, names_from = Rcombi)
    
write_csv(out_std_Rsum, paste0(tabledir, "/out_std_Rsum.csv"))
    
    
pt_phi_Ralleles_plt <- ggplot(caseH_phi, aes(x=CH1, y=CP1)) + 
  geom_point(aes(color=factor(R_alleles, levels = c(rev(as.character(0:3))))), show.legend=T, size = 0.3) +
  facet_grid(cols= vars(Species_label), rows = vars(phi), labeller=label_bquote(rows = phi == .(as.character(phi)))) +
  scale_color_manual(
    values = c("0" = "burlywood3", 
                "1" = "#FF8F00",
                "2" = "#64B5F6",
                "3" = "#3C5488FF"), 
    labels = c("0" = "none", 
                "1" = "private only", 
                "2" = "ancestral only", 
                "3" = "both"), 
    name = "Resistance alleles\nmaintained", 
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  theme_bw() +
  labs(x = bquote(bold("maximum cost of resistance"~Omega[H])),
       y = bquote(bold("maximum cost of virulence"~ Omega[P]))) +
  scale_x_continuous(limits = c(0, 0.3),
                     labels = c(0, 0.1, 0.2, 0.3),
                     expand = expansion(c(0,0), c(0,0.005)),
                     breaks = c(0, 0.1, 0.2, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3),
                     labels = c(0, 0.1, 0.2, 0.3),
                     expand = expansion(c(0,0), c(0,0.005)),
                     breaks = c(0, 0.1, 0.2, 0.3)) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 14),
        strip.background = element_rect(fill="white", color = "white"),
        strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")))  +
  theme(plot.margin = margin(0,0,0,0, 'cm'),
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.5, "cm"), 
        panel.spacing.y = unit(0.5, "cm")) 
    
f2 <- ggplot(caseH_phi , aes(x="", fill=factor(R_alleles, levels = c(rev(as.character(0:3)))))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
      y = "Proportion of 10,000 simulations") +
  theme_bw() +
  facet_grid(rows = vars(phi), cols = vars(Species_label), labeller=label_bquote(rows = phi == .(as.character(phi)) , cols = .(Species_label))) +
  scale_fill_manual(
    values = c("0" = "burlywood3", 
                "1" = "#FF8F00",
                "2" = "#64B5F6",
                "3" = "#3C5488FF"), 
    labels = c("0" = "none", 
                "1" = "private only", 
                "2" = "ancestral only", 
                "3" = "both"), 
    name = "", 
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(title = "", 
        x = "", 
        y = "Proportion out of\n10,000 simulations")+
  guides(fill="none") + 
  geom_text(data = Ralleles_summary , aes(y=prop_y, label = label, group = R_alleles), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 7,
            color = "white") +
  theme(plot.margin = margin(0,0,0,0, 'cm'))  + 
  theme(strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")), 
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
        panel.background = element_rect(fill = "gray60")) +
  theme(legend.text = element_text(size = 14))



# Summarize the R_alleles across all simulations
# This summary will be used to add the proportions to the generated barplot
# Note now label will be added if the proportion is < 5%
R_alleles_summary <- datH_dat %>% 
  group_by(phi_label, Species_label, R_alleles = factor(R_alleles, levels = rev(as.character(0:3)))) %>% 
  tally()  %>%
  group_by(phi_label, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  # If the proportion is smaller than 2.5% then make the label empty
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  )) |>
  mutate(Rcombi = case_when(
    R_alleles == 0 ~ "none",
    R_alleles == 1 ~ "private only",
    R_alleles == 2 ~ "ancestral only",
    R_alleles == 3 ~ "both", 
    TRUE ~ "unaccounted"
  ))

out_all_Rsum <- R_alleles_summary |>
  select(Species_label, phi_label, Rcombi, prop_y) |>
  mutate(prop_y = format(round(prop_y, digits =3), nsmall = 3)) |>
  pivot_wider(values_from = prop_y, names_from = Rcombi)

write_csv(out_all_Rsum, paste0(tabledir, "/out_all_Rsum.csv"))

# This is for the labelling the facets

# Generate the plot
R_alleles_acrossall_plt <- ggplot(datH_dat, aes(x="", fill= factor(R_alleles, levels = rev(as.character(0:3))))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(rows = vars(phi_label), cols = vars(Species_label), labeller=labeller(phi_label =label_parsed)) +
  scale_fill_manual(
    values = c("0" = "burlywood3", 
               "1" = "#FF8F00",
               "2" = "#64B5F6",
               "3" = "#3C5488FF"), 
    labels = c("0" = "none", 
               "1" = "private only", 
               "2" = "ancestral only", 
               "3" = "both"), 
    name = "Maintained\nresistance alleles", 
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(x = "", 
       y = "Proportion out of\n90,000 simulations") +
  theme_bw() +
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
  theme(legend.text = element_text(size = 13)) +
  #strip.background = element_rect(fill = "white", colour = "white")) +
  geom_text(data = R_alleles_summary, aes(y=prop_y, label = label), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 7,
            color = "white")


Fig2_new <- pt_phi_Ralleles_plt + theme(legend.position = "none") + ggtitle("Effect of maximum costs\nfor fixed shapes") + 
  plot_spacer() + 
  f2 + theme(legend.position = "none") +  ggtitle("Summary maximum\ncosts fixed shapes") +
  R_alleles_acrossall_plt +
  plot_layout(
    widths = c(2, 0.01, 1, 1),    # 1/3 : 2/3 column ratio
    heights = c(1)    # equal row heights
  ) +
  plot_annotation(
    tag_levels = list(c("(a)", "(b)", "(c)", "")),
    tag_prefix = "",
    tag_suffix = ""
  ) &
  theme(
    plot.tag = element_text(face = 2, size = 14),
    plot.title = element_text(face = 2, size = 16)
  )





pdf(paste0(main_figures,"/Fig_2_ralleles.pdf"), width = 15, height = 9)
print(Fig2_new)
dev.off()

    
# Figure 3

poly_summary <- caseH_phi %>% 
  group_by(phi, Species_label, poly = factor(poly, levels = rev(as.character(0:3)))) %>% 
  tally()  %>%
  group_by(phi, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))|>
  # If the proportion is smaller than 2.5% then make the label empty
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  )) |>
  mutate(poly = case_when(
    poly == 0 ~ "none",
    poly == 1 ~ "private only",
    poly == 2 ~ "ancestral only",
    poly == 3 ~ "both", 
    TRUE ~ "unaccounted"
  ))

out_std_poly <- poly_summary |>
  select(Species_label, phi, poly, prop_y) |>
  mutate(prop_y = format(round(prop_y, digits =3), nsmall = 3)) |>
  pivot_wider(values_from = prop_y, names_from = poly)

write_csv(out_std_poly, paste0(tabledir, "/out_std_poly.csv"))
    
    
    
pt_phi_poly_plt <- ggplot(caseH_phi, aes(x=CH1, y=CP1)) + 
  geom_point(aes(color=factor(poly, levels = rev(as.character(0:3)))), show.legend=T, size = 0.3) +
  facet_grid(cols= vars(Species_label), rows = vars(phi), labeller=label_bquote(rows = phi == .(as.character(phi)))) +
  scale_color_manual(values = c("0" = "burlywood3", 
                                "1" = "#FF8F00",
                                "2" = "#64B5F6",
                                "3" = "#3C5488FF"), 
                     labels = c("0" = "none", 
                                "1" = "private only", 
                                "2" = "ancestral only",
                                "3" = "both"), 
    name = "Maintained polymorphism", 
    drop = FALSE  # Ensures unused factor levels are retained
  )  +
  theme_bw() +
  labs(x = bquote(bold("maximum cost of resistance"~Omega[H])),
       y = bquote(bold("maximum cost of virulence"~ Omega[P]))) +
  scale_x_continuous(limits = c(0, 0.3), labels = c(0, 0.1, 0.2, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3), labels = c(0, 0.1, 0.2, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 14),
        strip.background = element_rect(fill="white", color = "white"),
        strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")))  +
  theme(plot.margin = margin(0,0,0,0, 'cm'),
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.5, "cm"), 
        panel.spacing.y = unit(0.5, "cm")) 
    
fpoly <- ggplot(caseH_phi , aes(x="", fill=factor(poly, levels = rev(as.character(0:3))))) +
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) +
  labs(x = "Species",
       y = "Proportion of 10,000 simulations") +
  theme_bw() +
  facet_grid(rows=vars(phi), cols = vars(Species_label),labeller=label_bquote(rows = phi == .(as.character(phi)) , cols = .(Species_label))) +
  scale_fill_manual(
    values = c("0" = "burlywood3",
               "1" = "#FF8F00",
               "2" = "#64B5F6",
               "3" = "#3C5488FF"),
    labels = c("0" = "none",
               "1" = "private only",
               "2" = "ancestral only",
               "3" = "both"), 
        name = "", 
    drop = FALSE) +
  labs(title = "",
       x = "",
       y = "Proportion out of\n10,000 simulations") +
  guides(fill="none") +
  geom_text(data = poly_summary , aes(y=prop_y, label = label, group = poly), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 6,
            color = "white") +
  theme(plot.margin = margin(0,0,0,0, 'cm'))  + 
  theme(strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")), 
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
        panel.background = element_rect(fill = "gray60")) +
  theme(legend.text = element_text(size = 14))

  
# Summarize the polymorphism data
# This is the summary to calcuate the proportions for each category
npoly_summary <- datH_dat %>% 
  group_by(phi_label, Species_label, poly) %>% 
  tally()  %>%
  group_by(phi_label, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))|>
  # If the proportion is smaller than 2.5% then make the label empty
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  )) |>
  mutate(poly = case_when(
    poly == 0 ~ "none",
    poly == 1 ~ "private only",
    poly == 2 ~ "ancestral only",
    poly == 3 ~ "both", 
    TRUE ~ "unaccounted"
  ))

out_all_poly <- npoly_summary |>
  select(Species_label, phi_label, poly, prop_y) |>
  mutate(prop_y = format(round(prop_y, digits =3), nsmall = 3)) |>
  pivot_wider(values_from = prop_y, names_from = poly)

write_csv(out_all_poly, paste0(tabledir, "/out_all_poly.csv"))


polyprops_acrossall_plt <- ggplot(datH_dat, aes(x="", fill= factor(poly, levels = rev(as.character(0:3))))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(rows = vars(phi_label), cols = vars(Species_label), labeller=labeller(phi_label =label_parsed)) +
  scale_fill_manual(
    values = c("0" = "burlywood3", 
               "1" = "#FF8F00",
               "2" = "#64B5F6",
               "3" = "#3C5488FF"), 
    labels = c("0" = "none", 
               "1" = "private only", 
               "2" = "ancestral only", 
               "3" = "both"), 
    name = "Maintained\npolymorphism", 
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(x = "", 
       y = "Proportion out of\n90,000 simulations") +
  theme_bw() +
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
  theme(legend.text = element_text(size = 13)) +
  geom_text(data = npoly_summary, aes(y=prop_y, label = label), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 6,
            color = "white") 



Fig3_new <- pt_phi_poly_plt + theme(legend.position = "none") + ggtitle("Effect of maximum costs\nfor fixed shapes") +
  plot_spacer() +
  fpoly + theme(legend.position = "none") +  ggtitle("Summary maximum\ncosts fixed shapes")  +
  polyprops_acrossall_plt +
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


pdf(paste0(main_figures,"/Fig_3_poly.pdf"), width = 15, height = 9)
print(Fig3_new)
dev.off()



# Figure 4 ----------------------------------------------------------------

scombined <- bind_rows(datH_dat, datP_dat)

## Summaries the ranges for the hosts across all simulations
nrange_summary <- scombined %>% 
  group_by(phi, Species, range_summary) %>% 
  tally()  %>%
  group_by(phi, Species) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(Species_label = case_when(
    Species == "H1" ~ "Host H",
    Species == "H2" ~ "Host M",
    Species == "P" ~ "Pathogen"
  ))|>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))




allrangebarplt <- ggplot(scombined |> mutate(Species_label = case_when(
  Species == "H1" ~ "Host H",
  Species == "H2" ~ "Host M",
  Species == "P" ~ "Pathogen"
)), aes(x="", fill=as.factor(range_summary))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(rows=vars(phi), cols = vars(Species_label), labeller=label_bquote(rows = phi == .(as.character(phi)) , cols = .(Species_label))) +
  scale_fill_manual(values = c("1" = "#AFB42B",
                               "2" = "#FFA726FF",
                               "3" = "#F57C00FF",
                               "4" = "#FF95A8FF",
                               "5" = "#EC407AFF",
                               "6" = "#C2185BFF",
                               "7" = "#8A4198FF",
                               "8" = "mediumorchid3",
                               "9" = "lightblue",
                               "10" = "darkseagreen2",
                               "11" = "limegreen",
                               "12" = "#008EA0FF",
                               "13" = "darkturquoise",
                               "14" = "#1A5355FF",
                               "15" = "black"),
                    labels = c("1" = "0",
                               "2" = "1",
                               "3" = "0 + 1",
                               "4" = "2",
                               "5" = "0 + 2",
                               "6" = "1 + 2",
                               "7" = "0 + 1 +2",
                               "8" = "3",
                               "9" = "0 + 3",
                               "10" = "1 + 3",
                               "11" = "0 + 1 +3",
                               "12" = "2 + 3",
                               "13" = "0 + 2 + 3",
                               "14" = "1 + 2 + 3",
                               "15" = "0 + 1 + 2  + 3"), drop = FALSE, name = "# of resistance (R) alleles\nor virulence (V) alleles/\n in maintained\nhaplotypes") +
  labs(title = "", 
       x = "", 
       y = "Proportion out of\n90,000 simulations") +
  theme_bw() +
  geom_text(data = nrange_summary , aes(y=prop_y, label = label, group = range_summary), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 6,
            color = "white")   +
  theme(strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")), 
        axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -0.75, size = 15),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm"), size = 15),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_text(face = "bold", size = 16),
        panel.spacing.x = unit(0, "cm"),
        panel.spacing.y = unit(0.5, "cm")) 

allrangebarplt <- tag_facet(allrangebarplt, x = 0.5, y = 1.05, size = 4, hjust = -0.5, vjust = 0) +
  theme(strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -0.75, size = 14),
        axis.title.y = element_text(size = 14),
        panel.grid = element_blank(),
        #strip.background = element_rect(fill = "gray90", colour = "gray40"),
        panel.background = element_rect(fill = "gray60"),
        panel.border = element_blank()) +
  guides(fill = guide_legend(ncol = 2))

# note default strip.background for theme_bw has fill = "gray85" and color = "gray20"

pdf(paste0(main_figures,"/Fig_4_ranges_across_all_shapes.pdf"), width = 8, height = 12)
print(allrangebarplt)
dev.off()

##  Supplementary figure


maintained_summary <- caseH_phi %>% 
  group_by(phi, Species_label, combi_maintained) %>% 
  tally()  %>%
  group_by(phi, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))

    
maintained_alleles_standard <- ggplot(caseH_phi, aes(x=CH1, y=CP1)) + 
      geom_point(aes(color= combi_maintained), show.legend=T, size = 0.05) + 
      facet_grid(cols = vars(Species_label), rows = vars(phi), labeller=label_bquote(rows = phi == .(as.character(phi)))) +
      theme_bw() + 
      scale_color_manual( values = c("gray40", "mediumvioletred","#925E9FFF","mediumseagreen","#FDAF91FF","#AD002AFF","cadetblue1","#0099B4FF","#00468BFF","lightyellow"),
                          name = "maintained alleles\nin population",drop =FALSE) +
      guides(color=guide_legend(override.aes = list(size = 5, shape = 15))) +
      labs(x = bquote(bold("maximum cost of resistance "*Omega[H])), 
           y = bquote(bold("maximum cost of virulence "*Omega[P]))) +
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
  scale_x_continuous(limits = c(0, 0.3),
                     labels = c(0, 0.1, 0.2, 0.3),
                     expand = expansion(c(0,0), c(0,0.005)),
                     breaks = c(0, 0.1, 0.2, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3),
                     labels = c(0, 0.1, 0.2, 0.3),
                     expand = expansion(c(0,0), c(0,0.005)),
                     breaks = c(0, 0.1, 0.2, 0.3)) 


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
  geom_text(data = maintained_summary , aes(y=prop_y, label = label, group = combi_maintained), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 5,
            color = "white") +
  theme(plot.margin = margin(0,0,0,0, 'cm'))  + 
  theme(strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")), 
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
        panel.background = element_rect(fill = "gray60")) +
  theme(legend.text = element_text(size = 14))


    

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

pdf(paste0(supplementary_figures,"/S4_alleles_maintained.pdf"), width = 15, height = 9)
print(SFig4_new)
dev.off()

# Now make a more detailed version of the plot. First generate the corresponding summary
# Note set the label threshold to 5% here
combi_summary_precise <- datH_dat %>% 
  mutate(phi = factor(phi)) |>
  group_by(CH2, CP2, phi, Species_label, combi_maintained) %>% 
  tally()  %>%
  group_by(CH2, CP2, phi, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))


# Generate the plot
Rcombi_costcombi_plt_hostM <- ggplot(datH_dat |> mutate(phi = factor(phi)) |> filter(Species_label == "Host M"), aes(x=phi, 
                                             fill= factor(combi_maintained))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Proportion of host H", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols = vars(CH2),rows = vars(CP2), labeller=label_bquote(cols = Xi[H] == .(CH2), rows = Xi[P] == .(CP2))) +
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
  #scale_x_discrete(labels = c("Host H.phi == 0.5" = "Host H", "Host M.phi == 0.5" = "Host M"))


Rcombi_costcombi_plt_hostH <- ggplot(datH_dat |> filter(Species_label == "Host H"), aes(x=as.factor(phi), 
                                                                                        fill= factor(combi_maintained))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Proportion of host H", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols = vars(CH2),rows = vars(CP2), labeller=label_bquote(cols = Xi[H] == .(CH2), rows = Xi[P] == .(CP2))) +
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
  ggtitle("Host species H")#+
#scale_x_discrete(labels = c("Host H.phi == 0.5" = "Host H", "Host M.phi == 0.5" = "Host M"))

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

pdf(paste0(supplementary_figures,"/SFig5_maintained_detail.pdf"),width = 15, height=8)
print(SFig5_new)
dev.off()


# Supplementary Figure 7 ----------------------------------------------------------------

## Check what type of polymorphism is maintained at the ancestral gene in both species
datH_trans <- datH_dat |>
  pivot_wider(id_cols = c(CH1, CH2, CP1, CP2, phi, seed), names_from = Species, values_from = poly, names_prefix = "poly_") |>
  mutate(transspecific = case_when(
    poly_H1%in%c(2,3) & poly_H2%in%c(2,3) ~ "trans-species",
    poly_H1%in%c(2,3) & !poly_H2%in%c(2,3) ~ "only in host H",
    !poly_H1%in%c(2,3) & poly_H2%in%c(2,3) ~ "only in host M",
    TRUE ~ "none"
  )) |>
  mutate(transspecific = factor(transspecific, levels = c("none","only in host M","only in host H", "trans-species")))

# Generate the summary
combi_summary <- datH_trans %>% 
  filter(phi<1) |>
  group_by(phi,CH2, CP2, transspecific = factor(transspecific)) %>% 
  tally()  %>%
  group_by(phi,CH2,CP2) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.025 ~ "",
    TRUE ~ label_raw
  ))



ancestral_poly <- ggplot(datH_trans |> filter(phi<1), aes(x=factor(phi), fill= factor(transspecific))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  facet_grid(cols = vars(CH2), rows = vars(CP2), labeller=label_bquote(cols = c[H]^{(2)} == .(CH2), rows = c[{P}]^{(2)} == .(CP2))) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  labs(x = bquote("Proportion of host H ("*phi*")"), 
       y = "Proportion") +
  scale_fill_manual(values = rev(c("#8E0152FF","#C51B7DFF","#F1B6DAFF","gray85")),name="Polymorphism at\nancestral R-gene") +
  theme_bw() +
  geom_text(data = combi_summary , aes(y=prop_y, label = label, group = transspecific), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 3,
            color = "black", fontface = 2) +
  theme(strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = 2),
        axis.text = element_text(size = 11),
        axis.title.x = element_text(vjust = -0.75, size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white")) 

pdf(paste0(main_figures,"/Fig_6_ancestral_poly.pdf"), width = 8, height = 5)
print(ancestral_poly)
dev.off()
  



