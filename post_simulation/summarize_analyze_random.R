require("tidyverse")
require("argparser")
require("ggpubr")
require("grid")
require("paletteer")
require("gridExtra")



args <- arg_parser("Script to analyze the maintained resistance alleles across costs")
args <- add_argument(args, arg = "--prefix", help ="Constant prefix for all simulations",
                     default = "randomcombs",
                     type = "character")
args <- add_argument(args, arg = "--suffix", 
                     help = "suffix for the simulations",
                     default = "2025-01-27")
args <- add_argument(args, arg = "--nreps", 
                     help = "number of simulations",
                     default = "10000",
                     type = "numeric")
args <- add_argument(args, arg = "--parent_dir", 
                     help = "directory which has all the results",
                     default = "results_random",
                     type = "character")
args <- add_argument(args, arg = "--cH_2", 
                     help = "value of CH2",
                     default = "3",
                     type = "numeric")
args <- add_argument(args, arg = "--cP_2", 
                     help = "value of CP2",
                     default = "3",
                     type = "numeric")
args <- add_argument(args, arg = "--rseed", 
                     help = "seed",
                     default = "1700",
                     type = "numeric")
args <- add_argument(args, arg = "--phi", 
                     help = "phi value",
                     default = "0.5",
                     type = "numeric")
args <- add_argument(args, arg = "--hinfo",
                     help = "host genotype info file",
                     default = "TSPolyR2/auxillaries/genotypesH_info.txt", 
                     type = "character")
args <- add_argument(args, arg = "--pinfo",
                     help = "host genotype info file",
                     default = "TSPolyR2/auxillaries/genotypesP_info.txt", 
                     type = "character")



pargs <- parse_args(args)
prefix <- pargs[["prefix"]]
suffix <- pargs[["suffix"]]
phi <- pargs[["phi"]]
cH_2 <- pargs[["cH_2"]]
cP_2 <- pargs[["cP_2"]]
parent_dir <- pargs[["parent_dir"]]
nreps <- pargs[["nreps"]]
seed_value <- pargs[["rseed"]]
hinfo_file <- pargs[["hinfo"]]
pinfo_file <- pargs[["pinfo"]]



rgenes_maintained <- function(sum_dat_host, haplocol){
  haplos <- sum_dat_host |> 
    select(!!haplocol) |> pull()
  a <- lapply(haplos, function(x){
    maintained_dat <- combmaintaineddat(x)
    found_ranges <- colSums(maintained_dat)
    
    possible_ranges <- 0:3
    
    range_boolean <- possible_ranges%in%found_ranges
    range_boolean
    
    range_weight <- 2^(possible_ranges)
    
    range_summary <- sum(range_boolean*range_weight)
    
    R_mat <- as.matrix(t(combmaintaineddat(x)))
    summary_Ralleles_present <- (colSums(R_mat) > 0)
    Poly_vec <- apply(R_mat,2, function(x){length(unique(x))-1})
    
    ret_dat <- tibble(range_summary = range_summary, R_1 = summary_Ralleles_present[1], 
                      R_2 = summary_Ralleles_present[2],
                      R_3 = summary_Ralleles_present[3],
                      poly_1 = Poly_vec[1],
                      poly_2 = Poly_vec[2],
                      poly_3 = Poly_vec[3])
    return(ret_dat)
  })
  daten <- bind_rows(a)
  combined <- bind_cols(sum_dat_host, daten) |>
    mutate(R_alleles = case_when(
      Species == "H1" ~ R_2 + 2*R_3,
      Species == "H2" ~ R_1 + 2*R_3,
      Species == "P" ~ R_1 + 2*R_2 + 4*R_3,
      TRUE ~ -999
    )
    ) |>
    mutate(poly = case_when(
      Species == "H1" ~ poly_2 + 2*poly_3,
      Species == "H2" ~ poly_1 + 2*poly_3, 
      Species == "P" ~ poly_1 + 2*poly_2 + 4*poly_3,
      TRUE ~ -999
    ))
  return(combined)
}


combmaintaineddat <- function(maintained_unparsed){
  maintained <- strsplit(maintained_unparsed, ", ")
  maintained_tabs <- sapply(maintained,function(x){
    # Split string into single loci
    split_haplos <- strsplit(x,"")
  })
  # Combine the list of maintained haplotypes for different cost combinations into a single data frame
  # Each column gives the allelic state of a single maintained haplotype
  # Each row is a single locus. Hence the number of of rows should correspond to the number of loci in the simulation.
  maintained_dat <- do.call(cbind.data.frame,maintained_tabs)
  # Convert all columns from string to numeric
  maintained_dat <- data.frame(apply(maintained_dat,2,as.numeric,drop=F))
  return(maintained_dat)
}

outdir <- paste0(parent_dir, "/", prefix, "_CH2_", cH_2, "_CP2_", cP_2, "_random_combis_", suffix, "_phi_", phi)
summarydir <- paste0(parent_dir, "/end_summaries")
enddir <- paste0(parent_dir, "/end_results")
combidir <- paste0(parent_dir, "/simulated_combinations")

if(!dir.exists(summarydir)){
  dir.create(summarydir, recursive = T)
}

if(!dir.exists(enddir)){
  dir.create(enddir, recursive = T)
}

if(!dir.exists(combidir)){
  dir.create(combidir, recursive = T)
}

fileprefix <- paste0("CH2_",cH_2,"_CP2_",cP_2,"_random_r_",nreps,"_RSEED_",seed_value)

daten <- read_delim(paste0(outdir,"/simulations/",fileprefix,"_end.txt"), 
                    delim = " ", col_names = F) 

#names(daten) <- c("origin","time",paste0("H1_",0:3), paste0("H2_",0:3), paste0("P_",0:7))
names(daten) <- c("time",paste0("H1_",0:3), paste0("H2_",0:3), paste0("P_",0:7), "origin")
write_tsv(daten, paste0(enddir,"/",fileprefix,"_phi_", phi, ".tsv"))

info <- read_tsv(paste0(outdir,"/simulations/",fileprefix,".tsv"))
names(info) <- c("CH1","CP1", "phi", "seed", "simulationID")

write_tsv(info, paste0(combidir,"/",fileprefix,"_phi_", phi, ".tsv"))

info <- info |>
  mutate(ID= sprintf("%05d",1:nrow(info)))

# Convert the dynamics to a tibble
daten <- as_tibble(daten) 
# daten <- daten |> 
#   mutate(origin = gsub(paste0("_",suffix,"_dynamics.txt"),"",gsub(paste0(prefix,"_"),"",basename(origin)))) |>
#   mutate(origin = gsub("cP_2","cP2", origin)) |>
#   mutate(origin = gsub("cH_2","cH2", origin)) |>
#   mutate(origin = gsub("cP_1","cP1", origin)) |>
#   mutate(origin = gsub("cH_1","cH1", origin)) 

#namen_mat <- str_split_fixed(daten[["origin"]],"_",18)

# Get the information on the genotypes for the hosts
hl_info <- read.table(hinfo_file, header=T,
                      colClasses=c("character","character","character",rep("numeric",times=10))) %>%
  as_tibble() %>%
  mutate(HostID=paste(Host,ID,sep="_"))

# Read the pathogen genotype info file
pl_info <- read.table(pinfo_file, header=T,
                      colClasses=c("character","character","character", rep("numeric",times=10))) %>%
  as_tibble() %>%
  mutate(ParaID=paste(Para,ID,sep="_"))

# Convert the tibble into long format
daten_long <- daten %>% 
  pivot_longer(!c(time, origin),names_to="Type",values_to="freq") %>% 
  separate(Type,c("Species","ID"),sep="_",remove=F)

# Add the host genotype information
daten_long <- left_join(daten_long, hl_info %>% select(HostID, Genotype),
                        by=c("Type"="HostID"))

# Add the pathogen genotype information
daten_long <- left_join(daten_long, pl_info %>% select(ParaID, Genotype),
                        by=c("Type"="ParaID"))

# Replace the missing values
daten_long <- daten_long %>% 
  replace_na(list(Genotype.x = "", Genotype.y=""))
# Combine the two separate genotype columns into a single column
daten_long <- daten_long %>% 
  unite("Genotype",Genotype.x,Genotype.y,sep="") |>
  mutate(origin = basename(origin))

#namen_mat <- str_split_fixed(daten_long[["origin"]],"_",18)

daten_long <- daten_long |>
  mutate(ID=gsub("simulation_","",origin))

daten_long <- daten_long |>
  left_join(info, by="ID")

daten_long <- daten_long |> filter(ID%in%sprintf("%05d",1:nrow(info)))

# daten_long <- daten_long |>
#   mutate(phi = as.numeric(namen_mat[, which(namen_mat[1,] == "phi")+1])) |>
#   mutate(CH2 = as.numeric(namen_mat[, which(namen_mat[1,] == "cH2")+1])) |>
#   mutate(CP2 = as.numeric(namen_mat[, which(namen_mat[1,] == "cP2")+1])) |> 
#   mutate(CH1 = as.numeric(namen_mat[, which(namen_mat[1,] == "cH1")+1])) |>
#   mutate(CP1 = as.numeric(namen_mat[, which(namen_mat[1,] == "cP1")+1])) |>
#   mutate(CH2 = as.numeric(namen_mat[, which(namen_mat[1,] == "cH2")+1])) |>
#   mutate(seed = as.numeric(namen_mat[, which(namen_mat[1,] == "seed")+1])) 



daten_long_maintained <- daten_long |>
  filter(freq > 0)

# daten_long_maintained <- daten_long_maintained |>
#   group_by(phi, CH2, CP2, CH1, CP1, seed, Species) |>
#   mutate(genotypes_remain = paste(Genotype, collapse=", ")) |>
#   ungroup()

daten_long_maintained <- daten_long_maintained |>
  group_by(phi, CH1, CP1, seed, Species) |>
  summarize(genotypes_remain = paste(Genotype, collapse=", ")) |>
  ungroup()

out <- rgenes_maintained(daten_long_maintained, "genotypes_remain")

outP <- out |> 
  filter(Species == "P") |>
  mutate(range_summary = factor(range_summary, levels = as.character(1:15)))

sizept <- 0.9

prangeplt <- ggplot(outP, aes(CH1, CP1)) + 
  geom_point(aes(color=range_summary), show.legend = T, size = sizept) +
  scale_color_manual(values = c("1" = "#AFB42B",
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
                                "15" = "0 + 1 + 2  + 3"), drop = FALSE, name = "# virulence alleles/\nhaplotype") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(x = bquote(c[H]^{(1)}),
       y = bquote(c[P]^{(1)})) 


pbarplt <- ggplot(outP, aes(x=Species, fill=range_summary)) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols=vars(phi), labeller=label_bquote(cols = phi == .(phi))) +
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
                                "15" = "0 + 1 + 2  + 3"), drop = FALSE, name = "# virulence alleles/\nhaplotype") +
  labs(title = "# virulence alleles in maintained haplotypes", 
       x = "Host species", 
       y = "Proportion") +
  theme_minimal() +
  labs(x = "Pathogen",
       y = "Proportion") 

write_tsv(outP |> mutate(CH2=cH_2, CP2=cP_2), paste0(summarydir,"/outP_", fileprefix , "_phi_", phi,".tsv"))

out_test <- out |>
  filter(Species%in%c("H1","H2")) |>
  mutate(R_alleles = factor(R_alleles, levels = c("0", "1", "2", "3"))) |>
  mutate(poly = factor(poly, levels = c("0", "1", "2", "3"))) |>
  mutate(range_summary = factor(range_summary, levels = c("1", "2", "3", "4","5", "6", "7"))) |>
  #filter(!(Species=="H2" & phi == 1) )
  filter(!(Species=="H2" & phi == 1) & !is.na(phi))

write_tsv(out_test |> mutate(CH2=cH_2, CP2=cP_2), paste0(summarydir,"/outH_", fileprefix , "_phi_", phi,".tsv"))

rbarplt <- ggplot(out_test, aes(x=Species, fill=R_alleles)) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols=vars(phi), labeller=label_bquote(cols = phi == .(phi))) +
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
  labs(title = "Maintained R alleles", 
       x = "Host species", 
       y = "Proportion") +
  theme_minimal() #+
  #scale_y_continuous(labels = scales::percent) +  # Formatting y-axis as percentage
  #theme(legend.position = "bottom")


polybarplt <- ggplot(out_test, aes(x=Species, fill=poly)) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  facet_grid(cols=vars(phi), labeller=label_bquote(cols = phi == .(phi))) +
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
  labs(title = "Maintained polymorphism", 
       x = "Host species", 
       y = "Proportion") +
  theme_minimal() 


rangebarplt <- ggplot(out_test, aes(x=Species, fill=range_summary)) + 
  geom_bar(stat = "count", position = "fill", show.legend = TRUE) + 
  labs(x = "Species", 
       y = "Proportion") +
  theme_minimal() +
  #facet_grid(cols=vars(phi), rows = vars(CP1), labeller=label_bquote(cols = phi == .(phi))) +
  facet_grid(cols=vars(phi), labeller=label_bquote(cols = phi == .(phi))) +
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
    name = "", 
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  labs(title = "# of NLRs/haplotype", 
       x = "Host species", 
       y = "Proportion") +
  theme_minimal() 


rangeplt <- ggplot(out_test, aes(x=CH1, y=CP1)) + 
  geom_point(aes(color=range_summary), show.legend=T, size = sizept ) +
  facet_grid(cols = vars(Species)) +
  scale_color_manual(
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
    name = "# of NLRs/\nhaplotype", 
    drop = FALSE  # Ensures unused factor levels are retained
  ) +
  theme_bw() +
  labs(x = bquote(c[H]^{(1)}),
       y = bquote(c[P]^{(1)})) +
  scale_x_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
  #theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Set fixed margin around plot area
  #      legend.key.width = unit(1, "cm"), # Set legend key width
  #      legend.key.height = unit(1, "cm"))  # Set legend key height
  #ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(CH2)*","~c[P]^{(2)}==.(CH2)))
  #ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==3*","~c[P]^{(2)}==3))




rplt <- ggplot(out_test, aes(x=CH1, y=CP1)) + 
  geom_point(aes(color=R_alleles), show.legend=T, size = sizept ) +
  facet_grid(cols = vars(Species)) +
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
  labs(x = bquote(c[H]^{(1)}),
       y = bquote(c[P]^{(1)})) +
  scale_x_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
  #theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Set fixed margin around plot area
  #     legend.key.width = unit(1, "cm"), # Set legend key width
  #      legend.key.height = unit(1, "cm"))  # Set legend key height
  #ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(CH2)*","~c[P]^{(2)}==.(CH2)))
  #ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==3*","~c[P]^{(2)}==3))

polyplt <- ggplot(out_test , aes(x=CH1, y=CP1)) + 
  geom_point(aes(color=poly), show.legend=T, size = sizept ) +
  facet_grid(cols = vars(Species)) +
  scale_color_manual(
    values = c("0" = "#B2DFDBFF", 
               "1" = "#00A087FF",
               "2" = "#64B5F6",
               "3" = "#3C5488FF"), 
    labels = c("0" = "none", 
               "1" = "private only", 
               "2" = "ancestral only", 
               "3" = "both"), 
    name = "maintained\npolymorphism", 
    drop = FALSE  # Ensures unused factor levels are retained
  )  +
  theme_bw() +
  labs(x = bquote(c[H]^{(1)}),
       y = bquote(c[P]^{(1)})) +
  scale_x_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  scale_y_continuous(limits = c(0, 0.3), expand = expansion(c(0,0), c(0,0.005)), breaks = c(0, 0.1, 0.2, 0.3)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
  #theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Set fixed margin around plot area
  #      legend.key.width = unit(1, "cm"), # Set legend key width
  #      legend.key.height = unit(1, "cm"))  # Set legend key height
  #ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(CH2)*","~c[P]^{(2)}==.(CH2)))
  #ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==3*","~c[P]^{(2)}==3))

#outplt <- ggarrange(rplt, polyplt, rangeplt, nrow = 3) 
#annotate_figure(outplt, top = text_grob(bquote("Shape parameters:"~c[H]^{(2)}==3*","~c[P]^{(2)}==3), hjust = 0, 
#                                        x = unit(45, "pt")))

                                        #color = "red", face = "bold", size = 14))

#grid.draw(rbind(ggplotGrob(rplt), ggplotGrob(polyplt), ggplotGrob(rangeplt), ggplotGrob(prangeplt), size="last"))
grid.arrange(rplt,  rangeplt, polyplt, prangeplt, nrow =2, ncol =2)

pdf(paste0(outdir,"/rplt_", fileprefix , "_phi_", phi,".pdf"), width = 10, height = 5.1)
rplt_out <- rplt +
  ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(cH_2)*","~c[P]^{(2)}==.(cP_2)*"  "~phi==.(phi)))
print(rplt_out)
dev.off()

pdf(paste0(outdir,"/rangeplt_", fileprefix , "_phi_", phi,".pdf"), width = 10, height = 5.1)
rangeplt_out <- rangeplt +
  ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(cH_2)*","~c[P]^{(2)}==.(cP_2)*"  "~phi==.(phi)))
print(rangeplt_out)
dev.off()

pdf(paste0(outdir,"/polyplt_", fileprefix , "_phi_", phi,".pdf"), width = 10, height = 5.1)
polyplt_out <- polyplt +
  ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(cH_2)*","~c[P]^{(2)}==.(cP_2)*"  "~phi==.(phi)))
print(polyplt_out)
dev.off()

pdf(paste0(outdir,"/prangeplt_", fileprefix , "_phi_", phi,".pdf"), width = 6.5, height = 5.1)
prangeplt_out <- prangeplt +
  ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(cH_2)*","~c[P]^{(2)}==.(cP_2)*"  "~phi==.(phi)))
print(prangeplt_out)
dev.off()



pdf(paste0(outdir,"/r_gene_maintained_", fileprefix , "_phi_", phi,".pdf"), height = 7.5, width = 9)
grid.draw(rbind(ggplotGrob(rplt), ggplotGrob(polyplt), size = "last"))
dev.off()

pdf(paste0(outdir,"/ranges_maintained_", fileprefix, "_phi_", phi,".pdf"), height = 8, width = 9)
#grid.draw(rbind(ggplotGrob(rangeplt), ggplotGrob(prangeplt), size = "last"))
grid.arrange(rangeplt, prangeplt, nrow=2,
             layout_matrix = rbind(c(1, 1, 1, 1, 1, 1),
                                   c(2, 2, 2, 2, NA, NA)))
dev.off()


pdf(paste0(outdir,"/bar_r_gene_maintained_",fileprefix, "_phi_", phi,".pdf"), height = 5, width = 9)
print(ggarrange(rbarplt, polybarplt, legend ="bottom", common.legend = TRUE))
dev.off()

pdf(paste0(outdir,"/bar_ranges_maintained_", fileprefix, "_phi_", phi,".pdf"), height = 5, width = 7)
print(ggarrange(rangebarplt, pbarplt))
dev.off()


pdf(paste0(outdir,"/rbarplt_", fileprefix , "_phi_", phi,".pdf"), height = 6, width = 4.2)
rbarplt_out <- rbarplt +
  ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(cH_2)*","~c[P]^{(2)}==.(cP_2)*"  "~phi==.(phi)))
print(rbarplt_out)
dev.off()


pdf(paste0(outdir,"/polybarplt_", fileprefix , "_phi_", phi,".pdf"), height = 6, width = 4.2)
polybarplt_out <- polybarplt +
  ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(cH_2)*","~c[P]^{(2)}==.(cP_2)*"  "~phi==.(phi)))
print(polybarplt_out)
dev.off()

pdf(paste0(outdir,"/rangebarplt_", fileprefix, "_phi_", phi,".pdf"), height = 6, width = 4.2)
rangebarplt_out <- rangebarplt +
  ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(cH_2)*","~c[P]^{(2)}==.(cP_2)*"  "~phi==.(phi)))
print(rangebarplt_out)
dev.off()


pdf(paste0(outdir,"/pbarplt_", fileprefix , "_phi_", phi,".pdf"), height = 6, width = 3.5)
pbarplt_out <- pbarplt +
  ggtitle(bquote("Shape parameters:"~c[H]^{(2)}==.(cH_2)*","~c[P]^{(2)}==.(cP_2)*"  "~phi==.(phi)))
print(pbarplt_out)
dev.off()
