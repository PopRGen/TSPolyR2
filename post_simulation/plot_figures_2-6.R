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
                     default = "Figures",
                     type="character")
pargs <- parse_args(args)

# Parse the command line arguments
indir <- pargs[["parentdir"]]
outdir <- pargs[["outdir"]]
suffix <- pargs[["suffix"]]

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
    R_alleles == 0 ~ "private S, ancestral S",
    R_alleles == 1 & poly == 0 ~ "private: R, ancestral: S",
    R_alleles == 1 & poly == 1 ~ "private: R/S, ancestral: S",
    R_alleles == 2 & poly == 0 ~ "private: S, ancestral: R",
    R_alleles == 2 & poly == 2 ~ "private: S, ancestral: R/S",
    R_alleles == 3 & poly == 0 ~ "private: R, ancestral: R",
    R_alleles == 3 & poly == 1 ~ "private: R/S, ancestral: R",
    R_alleles == 3 & poly == 2 ~ "private: R, ancestral: R/S",
    R_alleles == 3 & poly == 3 ~ "private: R/S, ancestral: R/S",
    TRUE ~ "forgotten"
  )) |>
  mutate(combi_maintained = factor(combi_maintained, levels = 
                                     c("private S, ancestral S",
                                       "private: R, ancestral: S",
                                       "private: R/S, ancestral: S",
                                       "private: S, ancestral: R",
                                       "private: R, ancestral: R",
                                       "private: R/S, ancestral: R",                                      
                                       "private: S, ancestral: R/S",                                       
                                       "private: R, ancestral: R/S",
                                       "private: R/S, ancestral: R/S")))

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
  theme_bw() +
  theme(legend.position = "none") + 
  ylab("Probability density") +
  theme(plot.margin = margin(0,0,0,0, 'cm')) +
  xlab(bquote(c[H]^{(1)})) #+ 
#+
#ggtitle(bquote(c[H]^{(2)} == 3~c[P]^{(2)}==3 ~ phi == 0.5))

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
  theme_bw() +
  theme(legend.position = "none") + 
  ylab("Probability density") +
  theme(plot.margin = margin(0,0,0,0, 'cm')) +
  xlab(bquote(c[P]^{(1)})) #+
#ggtitle(bquote(c[H]^{(2)} == 3~c[P]^{(2)}==3 ~ phi == 0.5))

pdf(paste0(main_figures,"/Fig_2a_std_case_rdensity_CH2_3_CP2_3_phi_0.5.pdf"), width = 7, height = 4.5)
print(rdens_std_plt)
dev.off()

pdf(paste0(main_figures,"/Fig_2c_std_case_rdensity_CH2_3_CP2_3_phi_0.5.pdf"), width = 7, height = 4.5)
print(rdensCP_std_plt)
dev.off()


rpnt_std_plt <- ggplot(standard_caseH,  aes(x = CH1, CP1, color = as.factor(maxR))) +
  geom_point(alpha=0.6, size = 0.8) +
  scale_color_manual(values= c("0" = "#9EBCDA",
                               "1" = "#8C6BB1",
                               "2" = "#810F7C",
                               "3" = "#4D004B"), name = "maximum number of\nR alleles/haplotype")  +
  theme_bw() +
  ylab(bquote("Cost of virulence"~c[P]^{(1)})) +
  xlab(bquote("Cost of resistance"~c[H]^{(1)})) +
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 15))) +
  theme( axis.text = element_text(size = 11),
         axis.title.x = element_text(vjust = -0.75, size = 12),
         axis.title.y = element_text(size = 12))#+

pdf(paste0(main_figures,"/Fig_2b_std_case_rpnt_CH2_3_CP2_3_phi_0.5.pdf"), width = 7, height = 6)
print(rpnt_std_plt)
dev.off()


composite_nb_alleles  <- (rdens_std_plt + scale_x_continuous(position = "top")) + plot_spacer() + rpnt_std_plt + (rdensCP_std_plt + 
                                                                                                                    scale_y_continuous(position = "right") + 
                                                                                                                    scale_x_continuous(position = "top") + rotate()) +
  plot_layout(widths = c(4, 1), heights = c(1,4))  & 
  plot_annotation(tag_levels = list(c("(a)","(b)","(c)",""))) &
  theme(plot.tag = element_text(face = 2))

pdf(paste0(main_figures,"/Fig_2_std_case_composite_CH2_3_CP2_3_phi_0.5.pdf"), width = 8.5, height = 8.5)
print(composite_nb_alleles )
dev.off()



# Figure 3 and corresponding additional -------------------------

for(CH2_value in c(-3,0,3)){
  ## Generate the R-overviews
  for(CP2_value in c(-3,0,3)){
    
    caseH_phi <- datH_dat |>
      filter(CH2 == CH2_value & CP2== CP2_value) |>
      mutate(phi = factor(phi, levels=c(0.5,0.7,0.9,1)))
    
    
    caseP_phi <- datP_dat |>
      filter(CH2== CH2_value & CP2 == CP2_value) |>
      mutate(phi = factor(phi, levels=c(0.5,0.7,0.9,1)))
    
    
    Ralleles_summary <- caseH_phi %>% 
      group_by(phi, Species_label, R_alleles = as.factor(R_alleles)) %>% 
      tally()  %>%
      group_by(phi, Species_label) |>
      mutate(prop_y = n/sum(n)) |>
      mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
      mutate(label = case_when(
        prop_y < 0.025 ~ "",
        TRUE ~ label_raw
      ))
    
    
    
    pt_phi_Ralleles_plt <- ggplot(caseH_phi, aes(x=CH1, y=CP1)) + 
      geom_point(aes(color=as.factor(R_alleles)), show.legend=T, size = 0.3) +
      facet_grid(rows= vars(Species_label), cols = vars(phi), labeller=label_bquote(cols = phi == .(as.character(phi)))) +
      scale_color_manual(
        values = c("0" = "#B2DFDBFF", 
                   "1" = "#00A087FF",
                   "2" = "#64B5F6",
                   "3" = "#3C5488FF"), 
        labels = c("0" = "none", 
                   "1" = "private only", 
                   "2" = "ancestral only", 
                   "3" = "both"), 
        name = "R-alleles\nmaintained", 
        drop = FALSE  # Ensures unused factor levels are retained
      ) +
      theme_bw() +
      labs(x = bquote('maximum cost of resistance'~c[H]^{(1)}),
           y = bquote('maximum cost of virulence'~c[P]^{(1)})) +
      scale_x_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
      scale_y_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
      guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
      theme(axis.title = element_text(size = 12),
            strip.background = element_rect(fill="white", color = "white"),
            strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2))  +
      theme(plot.margin = margin(0,0,0,0, 'cm'),
            panel.grid = element_blank()) 
    
    f2 <- ggplot(caseH_phi , aes(x="", fill=factor(R_alleles))) + 
      geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
      labs(x = "Species", 
           y = "Proportion") +
      theme_bw() +
      facet_grid(cols=vars(phi), rows = vars(Species_label),labeller=label_bquote(cols = phi == .(as.character(phi)))) +
      scale_fill_manual(
        values = c("0" = "#B2DFDBFF", 
                   "1" = "#00A087FF",
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
           y = "Proportion of simulations")+
      guides(fill="none") + 
      geom_text(data = Ralleles_summary , aes(y=prop_y, label = label, group = R_alleles), 
                position = position_stack(vjust=0.5),  # Adjust label position
                size = 3,
                color = "darkorange", fontface = 2) +
      theme(plot.margin = margin(0,0,0,0, 'cm'))  + 
      theme(axis.title = element_text(size = 12), 
            panel.grid = element_blank(),
            axis.text.x = element_blank(),                 # 
            axis.ticks.x = element_blank(), 
            panel.spacing = unit(0, "mm"),                       # 
            strip.background = element_blank(),
            axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
            strip.text.y = element_blank(),
            strip.text.x = element_text(size = 12, face =2)
      ) +
      
      theme(panel.border = element_blank()) 
    
    Fig3ab <- f2 + pt_phi_Ralleles_plt + plot_layout(widths = c(2,5), guides = "collect") & 
      plot_annotation(tag_levels = list(c("(a)","(b)","(c)",""))) &
      theme(plot.tag = element_text(face = 2))
    
    
    
    if(CH2_value == 3 & CP2_value == 3 ){
      pdf(paste0(main_figures,"/Fig_3ab_cH1_cP1_variation_ralleles.pdf"), width = 12.5, height = 4.5)
      print(Fig3ab)
      dev.off()
    }else{
      pdf(paste0(supplementary_figures,"/Supp_for_Fig3ab_CH2_",CH2_value,"_CP2_",CP2_value,"_ralleles.pdf"), width = 12.5, height = 4.5)
      print(Fig3ab)
      dev.off()
    }
    
    poly_summary <- caseH_phi %>% 
      group_by(phi, Species_label, poly = as.factor(poly)) %>% 
      tally()  %>%
      group_by(phi, Species_label) |>
      mutate(prop_y = n/sum(n)) |>
      mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
      mutate(label = case_when(
        prop_y < 0.025 ~ "",
        TRUE ~ label_raw
      ))
    
    
    
    pt_phi_poly_plt <- ggplot(caseH_phi, aes(x=CH1, y=CP1)) + 
      geom_point(aes(color=as.factor(poly)), show.legend=T, size = 0.3) +
      facet_grid(rows= vars(Species_label), cols = vars(phi), labeller=label_bquote(cols = phi == .(as.character(phi)))) +
      scale_color_manual(
        values = c("0" = "#B2DFDBFF", 
                   "1" = "#00A087FF",
                   "2" = "#64B5F6",
                   "3" = "#3C5488FF"), 
        labels = c("0" = "none", 
                   "1" = "private only", 
                   "2" = "ancestral only", 
                   "3" = "both"), 
        name = "Maintained poly-\nmorphism", 
        drop = FALSE  # Ensures unused factor levels are retained
      ) +
      theme_bw() +
      labs(x = bquote('maximum cost of resistance'~c[H]^{(1)}),
           y = bquote('maximum cost of virulence'~c[P]^{(1)})) +
      scale_x_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
      scale_y_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
      guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
      theme(axis.title = element_text(size = 12),
            strip.background = element_rect(fill="white", color = "white"),
            strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2))  +
      theme(plot.margin = margin(0,0,0,0, 'cm'),
            panel.grid = element_blank()) 
    
    fpoly <- ggplot(caseH_phi , aes(x="", fill=factor(poly))) + 
      geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
      labs(x = "Species", 
           y = "Proportion") +
      theme_bw() +
      facet_grid(cols=vars(phi), rows = vars(Species_label),labeller=label_bquote(cols = phi == .(as.character(phi)))) +
      scale_fill_manual(
        values = c("0" = "#B2DFDBFF", 
                   "1" = "#00A087FF",
                   "2" = "#64B5F6",
                   "3" = "#3C5488FF"), 
        labels = c("0" = "none", 
                   "1" = "private only", 
                   "2" = "ancestral only", 
                   "3" = "both"), 
        name = "", 
        drop = FALSE  # Ensures unused factor levels are retained
      ) +
      geom_text(data = poly_summary , aes(y=prop_y, label = label, group = poly), 
                position = position_stack(vjust=0.5),  # Adjust label position
                size = 3,
                color = "darkorange", fontface = 2) +
      labs(title = "", 
           x = "", 
           y = "Proportion of simulations")+
      guides(fill="none") +
      theme(plot.margin = margin(0,0,0,0, 'cm'))  + 
      theme(axis.title = element_text(size = 12), 
            panel.grid = element_blank(),
            axis.text.x = element_blank(),                 # 
            axis.ticks.x = element_blank(), 
            panel.spacing = unit(0, "mm"),                       # 
            strip.background = element_blank(),
            axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
            strip.text.y = element_blank(),
            strip.text.x = element_text(size = 12, face =2)
      ) +
      theme(panel.border = element_blank()) 
    
    Fig3cd <- fpoly + pt_phi_poly_plt + plot_layout(widths = c(2,5), guides = "collect") & 
      plot_annotation(tag_levels = list(c("(a)","(b)","(c)",""))) &
      theme(plot.tag = element_text(face = 2))
    
    
    if(CH2_value == 3 & CP2_value == 3 ){
      pdf(paste0(main_figures, "/Fig_3cd_phi_variation_poly.pdf"), width = 12.5, height = 4.5)
      print(Fig3cd)
      dev.off()
    }else{
      pdf(paste0(supplementary_figures,"/Supp_for_Fig_3cd_CH2_",CH2_value,"_CP2_",CP2_value,"_poly.pdf"), width = 12.5, height = 4.5)
      print(Fig3cd)
      dev.off()
    }
    
    Fig3_combined <- (f2 + ggtitle("Maintained R-alleles")  + 
                        theme(title = element_text(face =2))) + (pt_phi_Ralleles_plt + guides(color="none") + ggtitle("Maintained R-alleles")  + 
                                                                   theme(title = element_text(face =2), 
                                                                         strip.text.y = element_text(face =1))) +
      (fpoly + ggtitle("Maintained polymorphism") + theme(title = element_text(face =2))) + (pt_phi_poly_plt + ggtitle("Maintained polymorphism") + 
                                                                                               theme(legend.title=element_blank(), title = element_text(face =2), 
                                                                                                     strip.text.y = element_text(face =1))) + 
      plot_layout(widths = c(2,5), guides = "collect", heights = c(1,1), nrow =2) & 
      plot_annotation(tag_levels = list(c("(a)","(b)","(c)","(d)"))) &
      theme(plot.tag = element_text(face = 2), legend.text = element_text(size = 14))
    
    
    if(CH2_value == 3 & CP2_value == 3 ){
      pdf(paste0(main_figures,"/Fig_3_standard_case.pdf"),width = 12.5, height = 8.5)
      print(Fig3_combined)
      dev.off()
    }else{
      pdf(paste0(supplementary_figures,"/Supp_for_Fig_3_combined_CH2_",CH2_value,"_CP2_",CP2_value,".pdf"), width = 12.5, height = 4.5)
      print(Fig3_combined)
      dev.off()
    }
    
    
    maintained_alleles_standard <- ggplot(datH_dat |> filter(CH2 == 3 & CP2==3), aes(x=CH1, y=CP1)) + 
      geom_point(aes(color=as.factor(combi_maintained)), show.legend=T, size = 0.2) + 
      facet_grid(rows= vars(Species_label), cols = vars(phi), labeller=label_bquote(cols = phi == .(as.character(phi)))) +
      theme_bw() + 
      scale_color_manual( values = c("gray40", "mediumvioletred","#925E9FFF","mediumseagreen","#FDAF91FF","#AD002AFF","#ADB6B6FF","#00468BFF","#0099B4FF","lightyellow"),
                          name = "maintained alleles in population",drop =FALSE) +
      guides(color=guide_legend(override.aes = list(size = 5, shape = 15))) +
      labs(x = bquote("maximum cost of resistance "*c[H]^{(1)}), 
           y = bquote("maximum cost of virulence "*c[P]^{(1)})) +
      theme_bw() +
      theme(strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 12, face = 2),
            axis.text = element_text(size = 11),
            axis.title.x = element_text(vjust = -0.75, size = 12),
            axis.title.y = element_text(size = 12),
            panel.grid = element_blank(),
            strip.background = element_rect(fill = "white", colour = "white")) 
    
    
    pdf(paste0(supplementary_figures,"/Supp_maintained_alleles_CH2_",CH2_value,"_CP2_",CP2_value,".pdf"), width = 10, height = 4.5)
    print(maintained_alleles_standard)
    dev.off()
    
    
  }}



# Figure 4 ----------------------------------------------------------------


# Summarize the R_alleles across all simulations
# This summary will be used to add the proportions to the generated barplot
# Note now label will be added if the proportion is < 2.5%
R_alleles_summary <- datH_dat %>% 
  group_by(phi_label, Species_label, R_alleles) %>% 
  tally()  %>%
  group_by(phi_label, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  # If the proportion is smaller than 2.5% then make the label empty
  mutate(label = case_when(
    prop_y < 0.025 ~ "",
    TRUE ~ label_raw
  ))

# This is for the labelling the facets
R_alleles_facet_labels <- R_alleles_summary |>
  group_by(phi_label) |>
  summarize(n=n()) |>
  mutate(flabel = paste0("(",c("I","II","III","IV"),")"))

# Generate the plot
R_alleles_acrossall_plt <- ggplot(datH_dat, aes(x=factor(Species_label), fill= factor(R_alleles))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols = vars(phi_label), labeller=labeller(phi_label =label_parsed)) +
  scale_fill_manual(
    values = c("0" = "#B2DFDBFF", 
               "1" = "#00A087FF",
               "2" = "#64B5F6",
               "3" = "#3C5488FF"), 
    labels = c("0" = "none", 
               "1" = "private only", 
               "2" = "ancestral only", 
               "3" = "both"), 
    name = "Maintained R-alleles/\nsimulation", 
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(x = "Host species", 
       y = "Proportion") +
  theme_bw() +
  theme(strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        axis.text = element_text(size = 11),
        axis.title.x = element_text(vjust = -0.75, size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white")) +
  geom_text(data = R_alleles_summary, aes(y=prop_y, label = label, group = R_alleles), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 4,
            color = "white") +
  geom_text(data = R_alleles_facet_labels, aes(x=0.7, y=1.05, label=flabel, group=phi_label),
            vjust = 0.5,
            inherit.aes = FALSE)

pdf(paste0(main_figures,"/Fig_4a_r_alleles_across_shapes.pdf"), width = 8, height = 5)
print(R_alleles_acrossall_plt)
dev.off()

# Summarize the polymorphism data
# This is the summary to calcuate the proportions for each category
npoly_summary <- datH_dat %>% 
  group_by(phi_label, Species_label, poly) %>% 
  tally()  %>%
  group_by(phi_label, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.025 ~ "",
    TRUE ~ label_raw
  ))

npoly_facet_labels <- npoly_summary |>
  group_by(phi_label) |>
  summarize(n=n()) |>
  mutate(flabel = paste0("(",c("I","II","III","IV"),")"))

polyprops_acrossall_plt <- ggplot(datH_dat, aes(x=factor(Species_label), fill= factor(poly))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_bw() +
  #facet_grid(rows=vars(CH2), cols = vars(phi), labeller=label_bquote(cols = phi == .(phi), rows = c[H]^{(2)} == .(CH2))) +
  #facet_grid(rows=vars(shape_combi), cols = vars(phi), labeller=label_bquote(cols = phi == .(phi))) +
  facet_grid(cols = vars(phi_label), labeller=labeller(phi_label =label_parsed)) +
  scale_fill_manual(
    values = c("0" = "#B2DFDBFF", 
               "1" = "#00A087FF",
               "2" = "#64B5F6",
               "3" = "#3C5488FF"), 
    labels = c("0" = "none", 
               "1" = "private only", 
               "2" = "ancestral only", 
               "3" = "both"), 
    name = "Maintained poly-\nmorphism simulation", 
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(title = "", 
       x = "", 
       y = "Proportion") +
  theme_bw()  + 
  theme(strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        axis.text = element_text(size = 11),
        axis.title.x = element_text(vjust = -0.75, size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white")) +
  geom_text(data = npoly_summary, aes(y=prop_y, label = label, group = poly), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 4,
            color = "white") +
  geom_text(data = npoly_facet_labels, aes(x=0.7, y=1.05, label=flabel, group=phi_label),
            vjust = 0.5,
            inherit.aes = FALSE)

pdf(paste0(main_figures,"/Fig_4b_poly_across_shapes.pdf"), width = 8, height = 5)
print(polyprops_acrossall_plt )
dev.off()

Fig4_new <- (R_alleles_acrossall_plt + guides(fill = "none") + 
               theme(plot.margin = margin(c(8,0,0,0)),
                     axis.title.x = element_blank(),
                     plot.title = element_text(face = 2, margin=margin(b = -2, unit = "pt"))) + 
               ggtitle("Maintained R-alleles")) + 
  (polyprops_acrossall_plt + theme(plot.margin = margin(c(8,0,0,0)),
                                   legend.text = element_text(size = 13),
                                   legend.title = element_blank(),
                                   plot.title = element_text(face = 2, margin=margin(b = -2, unit = "pt"))) +
     ggtitle("Maintained polymorphism")) + plot_layout(nrow = 2, guides = "collect") + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = 2))

pdf(paste0(main_figures,"/Fig_4_across_shapes.pdf"), width = 8, height = 8)
print(Fig4_new)
dev.off()

# Summary for the maintained alleles
combi_summary <- datH_dat %>% 
  group_by(phi_label, Species_label, combi_maintained) %>% 
  tally()  %>%
  group_by(phi_label, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.025 ~ "",
    TRUE ~ label_raw
  ))


Rcombi_acrossall_plt <- ggplot(datH_dat, aes(x=factor(Species_label), fill= factor(combi_maintained))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE, color = "white", linewidth = 0.1) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols = vars(phi_label),labeller=labeller(phi_label =label_parsed)) +
  scale_fill_manual(
    values = c("gray40","#925E9FFF","violetred","mediumseagreen","#AD002AFF","#FDAF91FF","lightblue","#0099B4FF","#00468BFF","lightyellow"),
    name = "maintained alleles in population",
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(x = "Host species", 
       y = "Proportion") +
  theme_bw() +
  theme(strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = 2),
        axis.text = element_text(size = 11),
        axis.title.x = element_text(vjust = -0.75, size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white")) +
  geom_text(data = combi_summary, aes(y=prop_y, label = label, group = combi_maintained), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 4,
            color = "white")

pdf(paste0(supplementary_figures,"/Addititional_alleles_maintained_acrossshapes.pdf"), width = 9, height = 4.5)
print(Rcombi_acrossall_plt)
dev.off()

# Now make a more detailed version of the plot. First generate the corresponding summary
# Note set the label threshold to 5% here
combi_summary_precise <- datH_dat %>% 
  group_by(CH2, CP2, phi_label, Species_label, combi_maintained) %>% 
  tally()  %>%
  group_by(CH2, CP2, phi_label, Species_label) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))


# Generate the plot
Rcombi_costcombi_plt <- ggplot(datH_dat, aes(x=interaction(factor(Species_label), phi_label), 
                                             fill= factor(combi_maintained))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE, color = "white", linewidth = 0.1) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols = vars(CH2),rows = vars(CP2), labeller=label_bquote(cols = c[H]^{(2)} == .(CH2), rows = c[{P}]^{(2)} == .(CP2))) +
  scale_fill_manual(
    values = c("gray40","#925E9FFF","violetred","mediumseagreen","#AD002AFF","#FDAF91FF","lightblue","#0099B4FF","#00468BFF","lightyellow"),
    name = "maintained alleles in population",
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(x = "Host species", 
       y = "Proportion") +
  theme_bw() +
  geom_text(data = combi_summary_precise, aes(y=prop_y, label = label, group = combi_maintained), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 2.5,
            color = "white") +
  theme(strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        legend.text = element_text(size = 11),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_text(size = 12, face = 2),
        axis.text = element_text(size = 11),
        axis.title.x = element_text(vjust = -0.75, size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
  scale_x_discrete(labels=c(rep(paste("Host",c("H","M")), times = 3), "Host H"))#+
  #scale_x_discrete(labels = c("Host H.phi == 0.5" = "Host H", "Host M.phi == 0.5" = "Host M"))


pdf(paste0(supplementary_figures,"/alleles_maintained_eachcombi.pdf"), width = 9, height = 4.5)
print(Rcombi_costcombi_plt)
dev.off()



# Figure 6 ----------------------------------------------------------------

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
  


# Figure 5 ----------------------------------------------------------------

## Summaries the ranges for the hosts across all simulations
nHrange_summary <- datH_dat %>% 
  group_by(phi, Species, range_summary) %>% 
  tally()  %>%
  group_by(phi, Species) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(Species_label = case_when(
    Species == "H1" ~ "Host H",
    Species == "H2" ~ "Host M"
  ))|>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.025 ~ "",
    TRUE ~ label_raw
  ))


allrangebarplt <- ggplot(datH_dat |> mutate(Species_label = case_when(
  Species == "H1" ~ "Host H",
  Species == "H2" ~ "Host M"
)), aes(x=factor(phi), fill=as.factor(range_summary))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  #facet_grid(rows=vars(CH2), cols = vars(phi), labeller=label_bquote(cols = phi == .(phi), rows =c[H]^ {(2)} == .(CH2))) +
  #facet_grid(rows=vars(shape_combi), cols = vars(phi), labeller=label_bquote(cols = phi == .(phi))) +
  facet_grid(cols = vars(Species_label)) +
  scale_fill_manual(
    values = c("1" = "#AFB42B",
               "2" = "#FFA726FF",
               "3" = "#F57C00FF",
               "4" = "#FF95A8FF",
               "5" = "#EC407AFF",
               "6" = "#C2185BFF",
               "7" = "#8A4198FF"), 
    labels = c("1" = "0",
               "2" = "1",
               "3" = "0 + 1",
               "4" = "2",
               "5" = "0 + 2",
               "6" = "1 + 2",
               "7" = "0 + 1 + 2"), 
    name = "# of R-alleles in\nmaintained haplotypes", 
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(title = "", 
       x = bquote("Proportion of Host H ("*phi*")"), 
       y = "Proportion") +
  theme_bw() +
  geom_text(data = nHrange_summary , aes(y=prop_y, label = label, group = range_summary), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 4,
            color = "white")   +
  coord_cartesian(xlim=c(1,4))
allrangebarplt <- tag_facet(allrangebarplt, x = 0.5, y = 1.05, size = 4, hjust = -0.5, vjust = 0) +
  theme(strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        axis.text = element_text(size = 11),
        axis.title.x = element_text(vjust = -0.75, size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white")) +
  guides(fill = guide_legend(ncol = 2))

# note default strip.background for theme_bw has fill = "gray85" and color = "gray20"

pdf(paste0(main_figures,"/Fig_5ab_hranges_across_all_shapes.pdf"), width = 8, height = 4)
print(allrangebarplt)
dev.off()


## Summaries the ranges for the pathogen across all simulations
nPrange_summary <- datP_dat %>% 
  group_by(phi, Species, range_summary) %>% 
  tally()  %>%
  group_by(phi, Species) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(Species_label = "Pathogen") |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.025 ~ "",
    TRUE ~ label_raw
  ))

nPrange_facet_labels <- nPrange_summary  |>
  group_by(phi) |>
  summarize(n=n()) |>
  mutate(flabel = paste0("(",letters[1:n()],")"))

allrangePbarplt <- ggplot(datP_dat |> mutate(Species_label = "Pathogen"), aes(x=factor(phi), fill=as.factor(range_summary))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  #facet_grid(rows=vars(CH2), cols = vars(phi), labeller=label_bquote(cols = phi == .(phi), rows =c[H]^ {(2)} == .(CH2))) +
  #facet_grid(rows=vars(shape_combi), cols = vars(phi), labeller=label_bquote(cols = phi == .(phi))) +
  facet_grid(cols = vars(Species_label)) +
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
                               "15" = "0 + 1 + 2  + 3"), drop = FALSE, name = "# of V alleles/\nhaplotype") +
  labs(title = "", 
       x = bquote("Proportion of Host H ("*phi*")"), 
       y = "Proportion") +
  theme_bw() +
  geom_text(data = nPrange_summary , aes(y=prop_y, label = label, group = range_summary), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 4,
            color = "white") +
  theme(strip.text = element_text(size = 12), 
        axis.text = element_text(size = 11),
        axis.title.x = element_text(vjust = -0.75, size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid = element_blank()) +
  #geom_text(data = nPrange_facet_labels, aes(x=as.numeric(factor(phi))-0.5, y=1.05, label=flabel, group=phi),
  #         hjust = 0,
  #         inherit.aes = FALSE,
  #         fontface = 2, size = 4) +
  annotate("text", x=0.5, y= 1.05, label = "(c)", fontface = 2, size = 4, hjust = -0.5, vjust = 0) +
  theme(strip.background = element_rect(fill="white", color = "white"),
        strip.text = element_text(size = 13, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2)) +
  guides(fill = guide_legend(ncol = 2)) + 
  coord_cartesian(xlim=c(1,4))


pdf(paste0(main_figures,"/Fig_5c_pranges_across_all_shapes.pdf"), width = 8, height = 5)
print(allrangePbarplt)
dev.off()

# Combine into a single plot
together_out <- (allrangebarplt) / (allrangePbarplt + plot_spacer())+ plot_layout(heights = c(3,3), 
                                                                                  widths = c(1,1), 
                                                                                  guides = 'collect',
                                                                                  axes = 'collect') 

pdf(paste0(main_figures,"/Fig_5_ranges_across_all_shapes_allthree.pdf"), width = 10, height = 8.5)
print(together_out)
dev.off()
