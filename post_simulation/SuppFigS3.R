#require("grid") # for doing some grid work
#require("gridExtra") # for simple ggplot layout tasks
require("patchwork") # needed for better plot layouts
require("egg") # needed for tag_facet
require("ggExtra") # for the ggMarginal function

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
plotdir <- "Figures/supplementary_figures"



standard_caseH <- datH_dat |>
  filter(XiH == 3 & XiP == 3 & phi == 0.5) 
standard_caseH <- standard_caseH |>
  mutate(maxR = sapply(standard_caseH[["range_summary"]], function(x){max(which(intToBits(x)==1))-1}))

standard_caseP <- datP_dat |>
  filter(XiH ==3 & XiP ==3 & phi == 0.5) 

standard_caseP <- standard_caseP |>
  mutate(maxV = sapply(standard_caseP[["range_summary"]], function(x){max(which(intToBits(x)==1))-1}))



rdens_std_plt <- ggplot(standard_caseH,  aes(x = OmegaH, color = as.factor(maxR), fill = as.factor(maxR))) +
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


rdensCP_std_plt <- ggplot(standard_caseH,  aes(x = OmegaP, color = as.factor(maxR), fill = as.factor(maxR))) +
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
  xlab(bquote("Maximum cost of virulence"~Omega[P]))



rpnt_std_plt <- ggplot(standard_caseH,  aes(x = OmegaH, y = OmegaP, color = as.factor(maxR))) +
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
         axis.title.y = element_text(size = 16))

SFig3 <- ggMarginal(rpnt_std_plt + 
                      scale_color_manual(values= c("0" = "#9EBCDA",
                                                   "1" = "#8C6BB1",
                                                   "2" = "#810F7C",
                                                   "3" = "#4D004B"),
                                         name = "maximum number of R alleles/haplotype") +
                      theme(axis.title.x = element_text(size = 20),
                            axis.title.y = element_text(size = 20),
                            legend.title = element_text(size = 18),
                            legend.text = element_text(size = 18), axis.text = element_text(size = 20)), 
                    type="density", groupColour = TRUE, groupFill = TRUE) 


pdf(paste0(plotdir,"/S3Fig.pdf"), width = 7, height = 6.5)
print(SFig3)
dev.off()

