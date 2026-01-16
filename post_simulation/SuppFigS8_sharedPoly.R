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
supplementary_figures <- "Figures/supplementary_figures"


## Check what type of polymorphism is maintained at the ancestral gene in both species
datH_trans <- datH_dat |>
  pivot_wider(id_cols = c(OmegaH, XiH, OmegaP, XiP, phi, seed), names_from = Species, values_from = poly, names_prefix = "poly_") |>
  mutate(transspecific = case_when(
    poly_H1%in%c(2,3) & poly_H2%in%c(2,3) ~ "trans-species",
    poly_H1%in%c(2,3) & !poly_H2%in%c(2,3) ~ "only in host H",
    !poly_H1%in%c(2,3) & poly_H2%in%c(2,3) ~ "only in host M",
    TRUE ~ "none"
  )) |>
  mutate(transspecific = factor(transspecific, levels = c("none","only in host M","only in host H", "trans-species")))

# Generate the summary
combi_summary <- datH_trans |>
  filter(phi<1) |>
  group_by(phi, XiH, XiP, transspecific = factor(transspecific)) |> 
  tally()  |>
  ungroup() |>
  group_by(phi, XiH, XiP) |>
  mutate(prop_y = n/sum(n)) |>
  mutate(label_raw = paste0(round(prop_y, digits=2)*100,"%")) |>
  mutate(label = case_when(
    prop_y < 0.05 ~ "",
    TRUE ~ label_raw
  ))



ancestral_poly <- ggplot(datH_trans |> filter(phi<1), aes(x=factor(phi), fill= factor(transspecific))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  facet_grid(cols = vars(XiH), rows = vars(XiP), labeller=label_bquote(cols = Xi[H] == .(XiH), rows = Xi[{P}] == .(XiP))) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  labs(x = bquote("Proportion of host H ("*phi*")"), 
       y = "Proportion") +
  scale_fill_manual(values = rev(c("#8E0152FF","#C51B7DFF","#F1B6DAFF","gray85")),name="Polymorphism at\nancestral R-gene") +
  theme_bw() +
  geom_text(data = combi_summary , aes(y=prop_y, label = label, group = transspecific), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 4,
            color = "black", fontface = 2) +
  theme(strip.text = element_text(size = 14, hjust = 0, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt"), face  = 2), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = 2),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -0.75, size = 14),
        axis.title.y = element_text(size = 14),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white")) 

pdf(paste0(supplementary_figures,"/S7Fig_ancestral_poly.pdf"), width = 8, height = 5)
print(ancestral_poly)
dev.off()


ggsave(filename = paste0(supplementary_figures,"/S7Fig_ancestral_poly.png"), 
       plot = ancestral_poly,
       width = 8, 
       height = 5, 
       units = "in", dpi = 400)