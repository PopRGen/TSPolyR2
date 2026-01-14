require("ggpubr") # for using ggtext, gradient_fill, gradient_color etc. and to add lm results etc
require("grid") # for doing some grid work
require("paletteer") # for more colors
require("gridExtra") # for simple ggplot layout tasks
require("ggridges") # if ridges to be plotted
require("viridis") # needed for access to color blind friendly palettes
require("hrbrthemes") # needed for access to additional themes
require("patchwork") # needed for better plot layouts
require("egg") # needed for tag_facet
require("xtable") # needed for exporting a formatted latex table
require("RColorBrewer") # needed for color coding the cell of the formatted latex table

# Description: 
# This script will plot the proportion of simulations for which genotype with given numbers of resistance and or virulence
# alleles are maintained
# For example: 0 + 1 indicates that the haplotype wi


indir <- "results_random"

# Results for the hosts and the pathogen
datH_dat <- read_tsv(paste0(indir, "/datH_dat.tsv"))
datP_dat <- read_tsv(paste0(indir, "/datP_dat.tsv"))

# Set the output directory
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

scombined <- bind_rows(datH_dat, datP_dat)  |> 
  mutate(Species_label = case_when(
    Species == "H1" ~ "Host H",
    Species == "H2" ~ "Host M",
    Species == "P" ~ "Pathogen"
))

## Summaries the ranges for the hosts across all simulations
nrange_summary <- scombined |>
  group_by(phi, Species, range_summary) |>
  tally() |>
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




allrangebarplt <- ggplot(scombined, aes(x = "", fill = as.factor(range_summary))) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "", y = "Proportion") +
  theme_minimal() +
  facet_grid(rows=vars(phi), cols = vars(Species_label), 
             labeller=label_bquote(rows = phi == .(as.character(phi)) , cols = .(Species_label))) +
  # Generate the fill scale
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
                               "7" = "0 + 1 + 2",
                               "8" = "3",
                               "9" = "0 + 3",
                               "10" = "1 + 3",
                               "11" = "0 + 1 + 3",
                               "12" = "2 + 3",
                               "13" = "0 + 2 + 3",
                               "14" = "1 + 2 + 3",
                               "15" = "0 + 1 + 2  + 3"), drop = FALSE,
                    name = "# of resistance (R) alleles\nor virulence (V) alleles\nin maintained\ngenotypes") +
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

# Add the visual aids 
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

ggsave(filename = paste0(main_figures,"/Fig_4_ranges_across_all_shapes.png"), 
       plot = allrangebarplt,
       width = 8, 
       height = 12, 
       units = "in", dpi = 400)


####### Export all values for different combinations of shape parameters and relative proportions of hosts

# Generate detailed overview for supplementary information
geno_summary_shape <- scombined |>
  mutate(genocombi = case_when(
    range_summary == "1" ~ "0",
    range_summary == "2" ~ "1",
    range_summary == "3" ~ "0 + 1",
    range_summary == "4" ~ "2",
    range_summary == "5" ~ "0 + 2",
    range_summary == "6" ~ "1 + 2",
    range_summary == "7" ~ "0 + 1 +2",
    range_summary == "8" ~ "3",
    range_summary == "9" ~ "0 + 3",
    range_summary == "10" ~ "1 + 3",
    range_summary == "11" ~ "0 + 1 + 3",
    range_summary == "12" ~ "2 + 3",
    range_summary == "13" ~ "0 + 2 + 3",
    range_summary == "14" ~ "1 + 2 + 3",
    range_summary == "15" ~ "0 + 1 + 2  + 3",
    TRUE ~ "unaccounted")) |>
  group_by(phi, Species_label, shape_combi, range_summary = factor(range_summary, levels = as.character(1:15))) |>
  # Keep the information for Xih, XiP, the resistance allele combination
  # Note the code is purpously written this way, as the summarize function
  # will throw an error if there is not single unique value as expected
  summarize(n = n(),
            XiH = unique(XiH),
            XiP = unique(XiP),
            genocombi = unique(genocombi))   |>
  # Calculate the proportion of simulations falling into each category
  mutate(prop_y = round(n/sum(n), digits= 4)) |>
  ungroup() |>
  # Generate combined label to help with converting from long to wide data format 
  # and preventing two separate rows (one for each host species) for a single combination of XiH, XiP and phi
  mutate(specgeno = paste(Species_label, genocombi, sep = ": ")) 


geno_summary_shape_hH <- geno_summary_shape |>
  ungroup() |>
  filter(Species_label == "Host H") |>
  select(genocombi, phi, prop_y, shape_combi, XiH, XiP) |>
  pivot_wider(values_from = prop_y, names_from = genocombi, values_fill = 0) |>
  arrange(shape_combi, phi) |>
  select(-shape_combi) |>
  relocate(phi, .after = "XiP")

write_tsv(geno_summary_shape_hH, paste0(tabledir, "/SuppTab_ranges_hostH_combi.tsv"))

geno_summary_shape_hM <- geno_summary_shape |>
  ungroup() |>
  filter(Species_label == "Host M") |>
  select(genocombi, phi, prop_y, shape_combi, XiH, XiP) |>
  pivot_wider(values_from = prop_y, names_from = genocombi, values_fill = 0) |>
  arrange(shape_combi, phi) |>
  select(-shape_combi) |>
  relocate(phi, .after = "XiP")

write_tsv(geno_summary_shape_hH, paste0(tabledir, "/SuppTab_ranges_hostM_combi.tsv"))

geno_summary_shape_P <- geno_summary_shape |>
  ungroup() |>
  filter(Species_label == "Pathogen") |>
  select(genocombi, phi, prop_y, shape_combi, XiH, XiP) |>
  pivot_wider(values_from = prop_y, names_from = genocombi, values_fill = 0) |>
  arrange(shape_combi, phi) |>
  select(-shape_combi) |>
  relocate(phi, .after = "XiP")

write_tsv(geno_summary_shape_P, paste0(tabledir, "/SuppTab_ranges_hostM_combi.tsv"))


# Function to generate color codes based on values
get_cell_color <- function(value) {
  # Choose a color palette
  color_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))
  # Convert value to color (values between 0 and 1)
  farbe <- color_palette(100)[as.integer(value * 100) + 1]
  farbe_string <- paste(sprintf("%.3f",col2rgb(farbe) / 255), collapse = ",")
  color_string <- paste0("\\cellcolor[rgb]{", farbe_string, "}")
  complete_string <- paste(color_string, value)
  return(complete_string)
}



# Function to print the LaTeX table with color coding
ltable_generate <- function(geno_summary, cwidth = 1){
  
  # Assign the cell color 
  ltable <- geno_summary |>
    mutate(across(!c(XiP, XiH, phi), ~map_chr(., ~ get_cell_color(.)))) 
  # Convert to a matrix
  ltable_matrix <- as.matrix(ltable)
  
  # Start the table
  table_latex <- paste0("\\begin{table}[ht]\n\\centering\n\\begin{tabular}{|ccc|", paste0(rep(paste0("p{", cwidth, "cm}"), times = ncol(ltable_matrix)-3), collapse = ""), "|}\n\\hline\n")
  
  # Column names (headers)
  table_latex <- paste(table_latex, paste("\\textbf{", colnames(ltable_matrix), "}", collapse=" & "), "\\\\ \\hline\n", sep="")
  
  # Loop through each row of the matrix
  for (i in 1:nrow(ltable_matrix)) {
    row_latex <- paste(ltable_matrix[i,], collapse = " & ")
    # row_latex <- paste(row_latex, "\\\\ \\hline\n", sep="")
    row_latex <- paste(row_latex, "\\\\ \n", sep="")
    table_latex <- paste(table_latex, row_latex, sep="")
  }
  
  # Finish the table
  table_latex <- paste(table_latex, "\\hline\n\\end{tabular}\n\\caption{}\n\\end{table}")
  return(table_latex)
}

# Generate the LaTeX code
latex_codehH <- ltable_generate(geno_summary_shape_hH)
latex_codehM <- ltable_generate(geno_summary_shape_hM)
latex_codeP <- ltable_generate(geno_summary_shape_P, cwidth = 1)

# Print LaTeX code
writeLines(latex_codehH, paste0(tabledir,"/ranges_hostH_shapes.tex"))
writeLines(latex_codehM, paste0(tabledir,"/ranges_hostM_shapes.tex"))
writeLines(latex_codeP, paste0(tabledir,"/ranges_pathogen_shapes.tex"))



