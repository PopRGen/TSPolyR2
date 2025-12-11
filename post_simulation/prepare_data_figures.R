require("tidyverse")

# Description: Script to combine all end results across the 10,000 simulations each for different combinations
# of XiH, XiP and phi. For each of these combination the 10,000 simulations have been generated
# by randomly drawing OmegaH and OmegaP from the corresponding priors and randomly initializing genotype frequencies in
# each species


# Set the input directory
indir <- "results_random/"

setwd(indir)


# Create an overview of all files in end summaries
files_toread <- list.files("end_summaries", full.names =T)

# Read file for the end points for the host
to_readH <- files_toread[grep(glob2rx("outH*.tsv"), list.files("end_summaries"))]
# Read files for the end points for the pathogen
to_readP <- files_toread[grep(glob2rx("outP*.tsv"), list.files("end_summaries"))]

# Read all host results
datH <- lapply(to_readH, function(x){
  out <- read_tsv(x)
})

# Combine all host results
datH_dat <- bind_rows(datH)


# Add a couple of new columns to the host master file
datH_dat <- datH_dat |>
  # Note in the previus version of the submission XiH was named CH2, XiP was named CP2
  mutate(shape_combi=paste(CH2, CP2, sep="_")) |>
  mutate(combi_label = paste("atop(c[H]^{(2)} ==", CH2,",", "c[P]^{(2)} == ", CP2,")")) |>
  mutate(phi_label = paste("phi ==", phi)) |>
  mutate(Species_label = case_when(
    Species == "H1" ~ "Host H",
    Species == "H2" ~ "Host M",
    TRUE ~ "other"
  ))

datH_dat <- datH_dat |>
  # Summarize the maintained alleles
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
                                       "private: R/S, shared: R/S"))) |>
  # Changes the names of the parameters to the new notations
  rename(c("OmegaH" = "CH1",
           "OmegaP" = "CP1",
           "XiH" = "CH2",
           "XiP" = "CP2" ))


# Read the files for the pathogen
datP <- lapply(to_readP, function(x){
  out <- read_tsv(x)
})

# Combine all files for the pathogen
datP_dat <- bind_rows(datP)
# Add more information to the pathogen df
datP_dat <- datP_dat |>
  mutate(shape_combi=paste(CH2, CP2, sep="_")) |>
  mutate(Species_label="Pathogen") |>
  # To accomodate the changes of parameter names in the new version of the manuscript
  rename(c("OmegaH" = "CH1",
           "OmegaP" = "CP1",
           "XiH" = "CH2",
           "XiP" = "CP2" ))

# Check if everything looks right
# There should be a total of 630,000 lines in the combined host data frame
# Explanation: 
# 9 (XiH and XiP) combinations x 10,000 (random OmegaH, omegaP combinations) x 3 (phi values: 0.5, 0.7, 0.9) x 2 hosts = 540,000
# 9 (XiH and XiP) combinations x 10,000 (random OmegaH, omegaP combinations) x 1 (phi values:1 ) x 1 host = 90,000
# Total: 630,000
# There should be a total of 360,000 lines in the pathogen combined data frame
# 9 (cH2 and cP2) combinations x 10,000 (random cH1, cP1 combinations) x 4 (phi values: 0.5, 0.7, 0.9, 1.0) = 360,000
nrow(datH_dat)
nrow(datP_dat)

# Make sure that each combination has exactly 10,000 simulations
# Return values should be empty tibbles
datH_dat |> 
  group_by(shape_combi, phi, Species_label) |>
  summarize(n=n()) |> filter(n !=10000)

datP_dat |> 
  group_by(shape_combi, phi, Species_label) |>
  summarize(n=n()) |> filter(n !=10000)
###

# Output the data files for further processing to generate figures 2,3,4 and S3-S5 and S7.
# See the corresponding scripts for generating these figures
write_tsv(datH_dat, "datH_dat.tsv")
write_tsv(datP_dat, "datP_dat.tsv")