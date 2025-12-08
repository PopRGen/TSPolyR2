suppressPackageStartupMessages(require("tidyverse"))
suppressPackageStartupMessages(require("argparser"))
suppressPackageStartupMessages(require("paletteer"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("patchwork"))
suppressPackageStartupMessages(require("gridExtra"))


options(scipen = 999)

# Description: Script to generate the updated S6_Figure
# Note for backward compatibility and outlining the reasoning of some checks
# for the new version of the plot the plot to generate the original S6_Figure is still included

indir <- "Data_SuppFig6"

XiH_select <- 0
XiP_select <- 0

omegaH_select <- 0.2
omegaP_select <- 0.3
phi_select <- 0.5 

plotdir <- "Figures/supplementary_figures" 
prefix <- "ext_reps"
suffix <- "2025-01-07"


outprefix <- paste0(plotdir,"/CH1_", omegaH_select,"_CP1_", omegaP_select, "_CH2_", XiH_select, "_CP2_", XiP_select, "_phi_", phi_select)


# Find all simulation files across all seeds for the given combination of cH2 and cP2
dfiles <- list.files(indir)


# Extract the exact paramter combination based on the file name
param_combi <- lapply(dfiles, function(f){
  preparse <- gsub(paste0(prefix,"_|_", suffix, "|_long.tsv"),"",f)
  precursor <- gsub("cP_2", "CP2", preparse)
  precursor <- gsub("cP_1", "CP1", precursor)
  precursor <- gsub("cH_1", "CH1", precursor)
  precursor <- gsub("cH_2", "CH2", precursor)
  split_name <- strsplit(precursor,"_")[[1]]
  phi <- as.numeric(split_name[which(split_name == "phi") + 1])
  XiH <- as.numeric(split_name[which(split_name == "CH2") + 1])
  omegaH <- as.numeric(split_name[which(split_name == "CH1") + 1])
  XiP <- as.numeric(split_name[which(split_name == "CP2") + 1])
  omegaP <- as.numeric(split_name[which(split_name == "CP1") + 1])
  seed <- as.numeric(split_name[which(split_name == "seed") + 1])
  return(data.frame(phi=phi, 
                    omegaH=omegaH, 
                    omegaP=omegaP, 
                    XiH=XiH, 
                    XiP=XiP, 
                    seed=seed))
})

# Combine into a data.frame
param_combi_dat <- bind_rows(param_combi)

# Add a numeric ID column
param_combi_dat <- param_combi_dat |> 
  mutate(numericID=1:n())

# Check which simulations are for the given CH2, CP2, phi, CH1, CP1 combination to investigated
# and get their numeric ID
to_investigate <- param_combi_dat |>
  filter(phi == phi_select & omegaH == omegaH_select & omegaP == omegaP_select) |>
  pull(numericID)

# Read the dynamic files of interest
investigation_list <- lapply(paste0(indir,"/",dfiles[to_investigate]),
                             function(x){
                               daten <- read_tsv(x) |>
                                 mutate(Genotype = factor(Genotype, levels = c("000","100","010","001","110","101","011","111")))
                             })
# Add the respective parameter combination
for(i in 1:length(to_investigate)){
  investigation_list[[i]] <- bind_cols(investigation_list[[i]], param_combi[[to_investigate[i]]])
}

# First time intervals for which the averages and standard deviations for the different genotypes are two be calculated
firstint <- 20000
# Second time interval length for which the averages and standard deviations for the different genotypes are two be calculated
secint <- 1000

# Combine into a single data frame
investigation_dat <- bind_rows(investigation_list)

# Calculate in which time interval a given data point falls into
investigation_dat <- investigation_dat |>
  mutate(interval = paste0((time%/%firstint)*firstint,"-",((time%/%firstint)+1)*firstint),
         interval2 = paste0((time%/%secint)*secint,"-",((time%/%secint)+1)*secint))


interval_firstfac <- paste0(seq(0, 100000, firstint),"-",seq(firstint, 100000 + firstint, firstint))

# Make sure that these are encoded as factors with the right order for the plots
investigation_dat <- investigation_dat |>
  mutate(interval_level = factor(interval, levels = interval_firstfac)) |>
  mutate(interval2_level = factor(interval2, levels = paste0(seq(0, 100000, secint),"-", seq(secint, 100000 + secint, secint))))

investigation_dat <- investigation_dat |>
  mutate(broader_cat = case_when(
    Species%in%c("H1","H2") & Genotype%in%c("100","010") ~ "R_priv(100/010)",
    Species%in%c("H1","H2") & Genotype%in%c("101","011") ~ "R_double(101/011)",
    TRUE ~ Genotype
  )) |>
  mutate(broader_cat = factor(broader_cat, levels = c("000","100","010","001","R_priv(100/010)","101","011","110","R_double(101/011)","111")))


# The last time interval will be excluded later
# By design of the simulation the simulations all end at t=100,000
# Therefore, the last time interval only includes one single data point and thus mean and standard deviation are not meaningful there
lastfirst <- paste0("100000-",  100000 + firstint)
lastsecond <- paste0("100000-",  100000 + secint)

# Summarize for each seed, species, genotype combination
sum_investigation_dat <- investigation_dat |>
  group_by(interval_level, seed, Species, Genotype) |>
  summarize(mean_freq=mean(freq),
            sd_freq = sd(freq),
            seed = unique(seed),
            phi = unique(phi),
            omegaH = unique(omegaH),
            omegaP = unique(omegaP),
            n_obs=n())

testH1 <- expand.grid(interval_level = interval_firstfac, 
                      seed=investigation_dat |> pull(seed) |> unique(), 
                      Species="H1", 
                      Genotype = factor(c("000","010","001","011")))

testH2 <- expand.grid(interval_level = interval_firstfac, 
                      seed=investigation_dat |> pull(seed) |> unique(), 
                      Species="H2", 
                      Genotype = c("000","100","001","101"))

testP <- expand.grid(interval_level = interval_firstfac, 
            seed=investigation_dat |> pull(seed) |> unique(), 
            Species="P", 
            Genotype = c("000","100","010","001","110","101","011","111"))

test <- bind_rows(testH1, testH2, testP) |> mutate(Genotype = factor(Genotype, levels = c("000","100","010","001","110","101","011","111")))

experiment <- left_join(test, sum_investigation_dat, by=c("interval_level", "seed", "Species", "Genotype")) |> 
  arrange(seed, Species, Genotype, interval_level)

# This a sanity check. 
print("Running sanity check 1.")
no_dpts <- experiment |> filter(is.na(n_obs)) |> 
  group_by(interval_level, seed) |> 
  summarize(n=n())

print("Running sanity check 2.")
# Check that everything worked fine
no_dpts |> filter(n!=16)


dpts <- experiment |> filter(!is.na(n_obs))

# If there are intervals without data points
if(nrow(no_dpts)>0){
  for(i in 1:nrow(no_dpts)){
    
    interval_select <- no_dpts[["interval_level"]][i]
    seed_select <- no_dpts[["seed"]][i]
    n_obs_prev <- 0
    interval_select_temp <- interval_select
    n_obs_next <- 0
    interval_select_next_temp <- interval_select
    
    while(n_obs_prev < 1){
      interval_prev <- interval_firstfac[which(interval_firstfac == interval_select_temp) - 1]
      n_obs_prev <- sum_investigation_dat |> 
        filter(interval_level == interval_prev  & seed == seed_select) |>
        pull(n_obs) |> unique()
      if(length(n_obs_prev) == 0){ n_obs_prev <- 0}
      interval_select_temp <- interval_prev
    }
    
    investigation_dat_prev <- investigation_dat |> 
      filter(interval_level == interval_prev  & seed == seed_select)
    mtime_prev <- investigation_dat_prev |>
      pull(time) |> 
      max()
    freqs_prev_last <- investigation_dat_prev |>
      filter(time == mtime_prev) |>
      select(Type, freq)
  
    while(n_obs_next < 1){
      interval_next <- interval_firstfac[which(interval_firstfac == interval_select_next_temp) +1]
      n_obs_next <- sum_investigation_dat |> 
        filter(interval_level == interval_next  & seed == seed_select) |>
        pull(n_obs) |> unique()
      if(length(n_obs_next) == 0){ n_obs_next <- 0}
      interval_select_next_temp <- interval_next
    }
    
    investigation_dat_next <- investigation_dat |> 
      filter(interval_level == interval_next  & seed == seed_select)
    mtime_next <- investigation_dat_next |>
      pull(time) |> 
      min()
    freqs_next_first <- investigation_dat_next |>
      filter(time == mtime_next)  |>
      select(Type, Genotype, Species, freq)
    
    compared <- full_join(freqs_prev_last, freqs_next_first, by = "Type", suffix = c("",".y")) |> 
      mutate(freq_diff = freq.y - freq) 
    
    check <- !any(compared[["freq_diff"]] > 1e-6)
    if(check){
      add_sum <- compared |>
        mutate(interval_level = factor(interval_select, levels = interval_firstfac)) |>
        rename("mean_freq" = "freq") |>
        mutate(seed = seed_select, sd_freq = 0, phi = phi_select, omegaH = omegaH_select, omegaP = omegaP_select, n_obs = -999) |>
        select(interval_level, seed, Species, Genotype, mean_freq, sd_freq, phi, omegaH, omegaP, n_obs)
    }else{
      stop("Something is wrong")
    }
    
    dpts <- bind_rows(dpts, add_sum)
    rm(add_sum)
  }
}

# Generate a translation table for positioning the calculations for the mean and standard
# deviation across seeds for a given time interval correctly
trans_dat <- tibble(Genotype=factor(levels(sum_investigation_dat[["Genotype"]]), 
                                    levels=levels(sum_investigation_dat[["Genotype"]]))) |> 
                      mutate(position_center=1:n())

# Summarize for each species, genotype combination across all seeds
sumover_investigation_dat <- investigation_dat |>
  group_by(interval_level, Species, Genotype) |>
  summarize(mean_freq = mean(freq),
            sd_freq = sd(freq),
            phi = unique(phi),
            omegaH = unique(omegaH),
            omegaP = unique(omegaP),
            n_obs = n())


filter_investigation_dat <- dpts |> 
  filter(as.character(interval_level) != lastfirst & mean_freq > 0)

# For a plot without the standard deviations/mean across simulations added see commit: ab73049

# Create a more complicated plot, which also contains, the results when averaging across all 
# seeds for the given interval 
# Note: filter_investigation_dat need to be the first layer to be added
# otherwise we will run into issues with the Genotype group correctly being located on the x-axis
sd_plt <- ggplot(data = filter_investigation_dat, 
       aes(Genotype, mean_freq, group=seed)) +
  # Now plot the point for each single seed
  geom_point(data=filter_investigation_dat, aes(color=Genotype), position=position_dodge(width=0.75)) +
  # At the error bar for each seed
  geom_errorbar(data=filter_investigation_dat, 
                aes(ymin = mean_freq - sd_freq, ymax = mean_freq + sd_freq, color=Genotype), 
                position=position_dodge(width=0.75)) +
  facet_grid(rows = vars(Species), 
             cols = vars(interval_level)) +
  theme_bw() +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  ylab("Mean genotype frequency") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, size = 8)) +
  labs(title = bquote(Omega[H]^{(1)}~"=" ~.(omegaH_select) ~ Omega[P]^{(1)} ~"="~ .(omegaP_select))) 

print(sd_plt)

# Check how many of the simulations got fixed for 011
fixcheck <- filter_investigation_dat |>
  filter(Species == "H1" & interval_level == "80000-100000" & Genotype == "011") |>
  mutate(fixed_011 = case_when(
    !(sd_freq == 0 & mean_freq == 1) ~ "no",
    TRUE ~ "yes"
  )) 

fcheck_result <- fixcheck |>
  group_by(fixed_011) |>
  summarize(n= n()) |>
  mutate(prop = n/sum(n)) |>
  mutate(ntoselect = round(prop*5, digits = 0)) |> 
  ungroup()

fixcheck <- left_join(fixcheck, fcheck_result, by = "fixed_011")

seedtplt <- fixcheck |>
  group_by(fixed_011) |>
  mutate(id=row_number()) |>
  filter(id <= ntoselect) |>
  pull(seed)

newS6_dat <- investigation_dat |> 
  filter(seed%in%seedtplt) |>
  mutate(species_label = case_when(
    Species == "H1" ~ "Host H",
    Species == "H2" ~ "Host M",
    Species == "P" ~ "Pathogen"
  )) |>
  mutate(ltest = paste0("rep", as.numeric(as.factor(seed)))) 


newS6_1 <- ggplot(newS6_dat, aes(time,freq)) + 
  facet_grid(cols = vars(species_label), rows = vars(ltest)) + 
  geom_line(aes(color = Genotype), linewidth = 0.8) + 
  #geom_point(aes(color = Genotype), size = 0.2) +
  theme_bw(base_size = 14) + 
  guides(color = guide_legend(override.aes = list(size = 3, linewidth = 2) ) ) +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  scale_x_continuous(limits = c(0, 5000), labels = c(0, 2500, 5000), breaks = c(0, 2500, 5000)) +
  coord_cartesian(xlim=c(0, 5000)) +
  xlab("time") +
  ylab("Genotype frequency") +
  ggtitle("Initial dynamics") +
  theme(axis.text = element_text(size = 14),
        axis.title.x = element_text(vjust = -0.75, size = 15),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm"), size = 15),
        panel.grid = element_blank(),
        title = element_text(face = "bold", size = 16),
        panel.spacing.y = unit(0.8, "cm"),
        panel.spacing.x = unit(0.8, "cm")) 

newS6_2 <- ggplot(newS6_dat, aes(time,freq)) + 
  facet_grid(cols = vars(species_label), rows = vars(ltest)) + 
  geom_line(aes(color = Genotype), linewidth = 0.8) + 
  #geom_point(aes(color = Genotype), size = 0.2) +
  theme_bw(base_size = 14) + 
  guides(color = guide_legend(override.aes = list(size = 3, linewidth = 2) ) ) +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  scale_x_continuous(limits = c(0, 100000), labels = c(0, 5, 100000), breaks = c(0, 50000, 100000)) +
  coord_cartesian(xlim=c(0,100000)) +
  xlab("time") +
  ylab("Genotype frequency") +
  ggtitle("Full length dynamics") +
  theme(axis.text = element_text(size = 14),
        axis.title.x = element_text(vjust = -0.75, size = 15),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm"), size = 15),
        panel.grid = element_blank(),
        title = element_text(face = "bold", size = 16),
        panel.spacing.y = unit(0.8, "cm"),
        panel.spacing.x = unit(0.8, "cm")) 



S6_new <- newS6_1 + newS6_2 +
  plot_layout(ncol = 2, guides = "collect") +
  # Add labels to the subfigures
  plot_annotation(
    tag_levels = list(c("(a)", "(b)")),
    tag_prefix = "",
    tag_suffix = ""
  ) &
  theme(legend.position = "bottom") &
  # Define the size of the subfigure labels and the plot subplot titles
  theme(
    plot.tag = element_text(face = 2, size = 14),
    plot.title = element_text(face = 2, size = 16)
  )

pdf("Figures/supplementary_figures/FigS6.pdf", width = 14, height = 8)
print(S6_new)
dev.off()