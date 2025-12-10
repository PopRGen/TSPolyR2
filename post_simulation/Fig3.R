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
require("ggExtra") # 
require("xtable") # needed for exporting a formatted latex table
require("RColorBrewer") # needed for color coding the cell of the formatted latex table

# Description: The script will generate a figure with three subfigures. 
# (a): The resistance gene polymorphism maintained for XiH=XiM=Xi_P=3 for different proportions of host H (rows) in each host species (columns).
# for the 10,000 randomly drawn combinations of Omega_H (x) and Omega_P 
# (b): The results from (a) aggregated as barplot for each phi (rows) x host species combination
# (c): The aggregated results across 9 combinations of shape parameters ($Xi_H = Xi_M \in \{-3, 0, 3\}$) x ($Xi_P \in \{-3, 0, 3\}$) and four different values of phi with 10,000 simulations for each combination.


XiH_value <- 3
XiP_value <- 3

indir <- "results_random"

datH_dat <- read_tsv(paste0(indir, "/datH_dat.tsv"))
datP_dat <- read_tsv(paste0(indir, "/datH_dat.tsv"))


outdir <- "Figures"

main_figures <- paste0(outdir,"/main_figures")
supplementary_figures <- paste0(outdir, "/supplementary_figures")
tabledir <- "upplementary_tables"

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

#########################################################
# Extract all results for the Xi_H = Xi_M = Xi_P = 3 
#########################################################


caseH_phi <- datH_dat |>
  filter(XiH == XiH_value & XiP == XiP_value) |>
  mutate(phi = factor(phi, levels=c(0.5, 0.7, 0.9, 1)))

caseP_phi <- datP_dat |>
  filter(XiH == XiH_value & XiP == XiP_value) |>
  mutate(phi = factor(phi, levels=c(0.5, 0.7, 0.9, 1)))

# Translations for R gene polymorphism are:
# 0: no polymorphism maintained
# 1: only private resistance gene polymorphism maintained and polymorphism at shared resistance gene lost
# 2: only shared resistance gene polymorphism maintained and polymorphism at private resistance gene lost
# 3: both resistance gene polymorphism are maintained

# Summarize the resistance alleles maintained for each host species and each value of phi
poly_summary <- caseH_phi |> 
  group_by(phi, Species_label, poly) |>
  tally()  |>
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
  )) |>
  mutate(poly = factor(poly, levels = c("none", "private only", "ancestral only", "both")))

# Output the essential parts of the summary table
out_std_poly <- poly_summary |>
  select(Species_label, phi, poly, prop_y) |>
  mutate(prop_y = format(round(prop_y, digits =3), nsmall = 3)) |>
  pivot_wider(values_from = prop_y, names_from = poly)

write_csv(out_std_poly, paste0(tabledir, "/out_std_poly.csv"))

# Generate subfigure (a)
pt_phi_poly_plt <- ggplot(caseH_phi, aes(x=OmegaH, y=OmegaP)) + 
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
  # Use theme_bw
  theme_bw() +
  # Define the axis labels
  labs(x = bquote(bold("maximum cost of resistance"~Omega[H])),
       y = bquote(bold("maximum cost of virulence"~ Omega[P]))) +
  # Define the axis-ticks 
  scale_x_continuous(limits = c(0, 0.3), labels = c(0, 0.1, 0.2, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3), labels = c(0, 0.1, 0.2, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  # Override the aesthetics for the color legend
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  # Refine the plot appearance
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 14),
        strip.background = element_rect(fill="white", color = "white"),
        strip.text = element_text(size = 16, hjust = 0.5, margin = margin(t=4.4, r=4.4, b=4.4, l=4.4, "pt")),
        plot.margin = margin(0,0,0,0, 'cm'),
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.5, "cm"), 
        panel.spacing.y = unit(0.5, "cm")) 

###################################
# Subfigure (b)
###################################

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
  # Remove the fill legend
  guides(fill="none") +
  geom_text(data = poly_summary , aes(x="", y=prop_y, label = label), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 6,
            color = "white",
            inherit.aes = FALSE) +
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


############################################################################
# Subfigure (c)
# Summarize the maintained polymorphism across all simulations from the different shape combinations
# This summary will be used to add the proportions to the generated barplot
# Note no label will be added if the proportion is < 5%
# ##########################################################################


# Summarize the polymorphism data
# This is the summary to calcuate the proportions for each category
npoly_summary <- datH_dat |>
  group_by(phi_label, Species_label, poly) |>
  tally()  |>
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
  ggtitle("Across cost function shapes\nand maximum costs") +
  theme(legend.text = element_text(size = 13)) +
  geom_text(data = npoly_summary, aes(y=prop_y, label = label), 
            position = position_stack(vjust=0.5),  # Adjust label position
            size = 6,
            color = "white") 

# Combine all three subfigures into a single figure

Fig3_new <- pt_phi_poly_plt + theme(legend.position = "none") + ggtitle("Effect of maximum costs\nfor fixed shapes") +
  plot_spacer() +
  fpoly + theme(legend.position = "none") +  ggtitle("Summary maximum\ncosts fixed shapes")  +
  polyprops_acrossall_plt +
  # Define the plot layout
  plot_layout(
    widths = c(2, 0.01,  1, 1),    # 1/3 : 2/3 column ratio
    heights = c(1)    # equal row heights
  ) +
  # Add labels to the subfigures
  plot_annotation(
    tag_levels = list(c("(a)", "(b)", "(c)")),
    tag_prefix = "",
    tag_suffix = ""
  ) &
  # Define the size of the subfigure labels and the plot subplot titles
  theme(
    plot.tag = element_text(face = 2, size = 14),
    plot.title = element_text(face = 2, size = 16)
  )

# Save the results to pdf
pdf(paste0(main_figures,"/Fig_3_poly.pdf"), width = 15, height = 9)
print(Fig3_new)
dev.off()


####### Export all values for different combinations of shape parameters and relative proportions of hosts

# Generate detailed overview for supplementary information
poly_summary_shape <- datH_dat |>
  mutate(polycombi = case_when(
    poly == 0 ~ "none",
    poly== 1 ~ "private only",
    poly == 2 ~ "ancestral only",
    poly == 3 ~ "both", 
    TRUE ~ "unaccounted"
  )) |>
  group_by(phi_label, Species_label, shape_combi, poly = factor(poly, levels = rev(as.character(0:3)))) |>
  # Keep the information for Xih, XiP, the resistance allele combination
  # Note the code is purpously written this way, as the summarize function
  # will throw an error if there is not single unique value as expected
  summarize(n = n(),
            XiH = unique(XiH),
            XiP = unique(XiP),
            polycombi = unique(polycombi))   |>
  # Calculate the proportion of simulations falling into each category
  mutate(prop_y = round(n/sum(n), digits= 4)) |>
  ungroup() |>
  # Generate combined label to help with converting from long to wide data format 
  # and preventing two separate rows (one for each host species) for a single combination of XiH, XiP and phi
  mutate(specpoly = paste(Species_label, polycombi, sep = ": ")) 


out_all_poly_shape <- poly_summary_shape |>
  ungroup() |>
  select(specpoly, phi_label, prop_y, shape_combi, XiH, XiP) |>
  pivot_wider(values_from = prop_y, names_from = specpoly, values_fill = 0) |>
  mutate(phi = gsub("phi == ", "", phi_label)) |>
  select( phi, XiH, XiP, starts_with("Host"), shape_combi) |>
  arrange(shape_combi, phi) |>
  select(-shape_combi) |>
  select(XiH, XiP, phi, `Host H: both`, `Host H: private only`, `Host H: ancestral only`, `Host H: none`, `Host M: both`, `Host M: private only`, `Host M: ancestral only`, `Host M: none`)

write_tsv(out_all_poly_shape, paste0(tabledir, "/SuppTab_poly_per_shapephi_combi.tsv"))



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
ltable_generate <- function(shape_summary){
  
  # Assign the cell color 
  ltable <- shape_summary |>
    mutate(across(starts_with("Host"), ~map_chr(., ~ get_cell_color(.)))) 
  # Convert to a matrix
  ltable_matrix <- as.matrix(ltable)
  
  
  #table_latex <- paste0("\\begin{table}[ht]\n\\centering\n\\begin{tabular}{", paste0(paste0("|",rep("p{1.5cm}", times = ncol(ltable_matrix))), collapse = ""), "|}\n\\hline\n")
  
  # Start the table
  table_latex <- paste0("\\begin{table}[ht]\n\\centering\n\\begin{tabular}{|ccc|", paste0(rep("p{1.5cm}", times = ncol(ltable_matrix)-3), collapse = ""), "|}\n\\hline\n")
  
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
  table_latex <- paste(table_latex, "\\end{tabular}\n\\caption{}\n\\end{table}")
  return(table_latex)
}

# Generate the LaTeX code
latex_code <- ltable_generate(out_all_poly_shape)

# Print LaTeX code
writeLines(latex_code, paste0(tabledir,"/poly_shapes.tex"))


