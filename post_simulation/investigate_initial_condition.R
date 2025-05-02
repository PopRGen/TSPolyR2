suppressPackageStartupMessages(require("tidyverse"))
suppressPackageStartupMessages(require("argparser"))
suppressPackageStartupMessages(require("paletteer"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("patchwork"))
suppressPackageStartupMessages(require("gridExtra"))


options(scipen = 999)


args <- arg_parser("Command line argument parser")
args <- add_argument(args,arg="--prefix",
                     help="Prefix to append to all simulations",
                     default = "testrun_small",
                     type="character")
args <- add_argument(args,arg="--suffix",
                     help="Suffix to append to all filenames",
                     type="character",
                     default="2024-11-15")
args <- add_argument(args,arg="--cH_2",
                     help="Value of cH_2 to use for the simulations",
                     default =3,
                     type="numeric")
args <- add_argument(args,arg="--cP_2",
                     help="Value of cP_2 to use for simulations",
                     default =0,
                     type="numeric")
args <- add_argument(args,arg="--indir",
                     help="Top level directory for the run",
                     default ="results_r50/all_dynamics_long",
                     type="character")
args <- add_argument(args, arg="--plotdir",
                     help="Directory for storing the simulations",
                     default="Figures_seed_comparison/extended_analyses",
                     type="character")
args <- add_argument(args, arg="--tselect",
                     help="Time from which onwards frequencies should be averaged",
                     default = 40000,
                     type="numeric")
args <- add_argument(args, arg="--phi",
                     help="Value of phi to analyze",
                     default = 0.5,
                     type="numeric")
args <- add_argument(args, arg="--cH_1",
                     help="Value of cH_1 to analyze",
                     default = 0.1,
                     type="numeric")
args <- add_argument(args, arg="--cP_1",
                     help="Value of cP_1 to analyze",
                     default = 0.2,
                     type="numeric")
args <- add_argument(args, "--nrep", 
                     help="number of different initial seed", 
                     default = 50,
                     type = "numeric")
args <- add_argument(args, "--seed_start", 
                     help="lower seed to start with", 
                     default = 400,
                     type = "numeric")
args <- add_argument(args, "--seed_interval", 
                     help="seed interval tp ise", 
                     default = 400,
                     type = "numeric")

# Parse command line arguments
pargs <- parse_args(args)
prefix <- pargs[["prefix"]]
suffix <- pargs[["suffix"]]
CH2 <- as.numeric(pargs[["cH_2"]])
CP2 <- as.numeric(pargs[["cP_2"]])
CH1_select <- as.numeric(pargs[["cH_1"]])
CP1_select <- as.numeric(pargs[["cP_1"]])
phi_select <- as.numeric(pargs[["phi"]])
indir <- pargs[["indir"]]
plotdir <- pargs[["plotdir"]]
tselect <- as.numeric(pargs[["tselect"]])
nrep <- pargs[["nrep"]]
seed_start <- pargs[["seed_start"]]
seed_interval <- pargs[["seed_interval"]]


if(!dir.exists(plotdir)){
  dir.create(plotdir)
}

outprefix <- paste0(plotdir,"/CH1_", CH1_select,"_CP1_", CP1_select, "_CH2_", CH2, "_CP2_", CP2, "_phi_", phi_select)


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
  CH2 <- as.numeric(split_name[which(split_name == "CH2") + 1])
  CH1 <- as.numeric(split_name[which(split_name == "CH1") + 1])
  CP2 <- as.numeric(split_name[which(split_name == "CP2") + 1])
  CP1 <- as.numeric(split_name[which(split_name == "CP1") + 1])
  seed <- as.numeric(split_name[which(split_name == "seed") + 1])
  return(data.frame(phi=phi, 
                    CH1=CH1, 
                    CP1=CP1, 
                    CH2=CH2, 
                    CP2=CP2, 
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
  filter(phi == phi_select & CP1 == CP1_select & CH1 == CH1_select) |>
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
            CH1 = unique(CH1),
            CP1 = unique(CP1),
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

no_dpts <- experiment |> filter(is.na(n_obs)) |> 
  group_by(interval_level, seed) |> 
  summarize(n=n())

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
        mutate(seed = seed_select, sd_freq = 0, phi = phi_select, CH1 = CH1_select, CP1 = CP1_select, n_obs = -999) |>
        select(interval_level, seed, Species, Genotype, mean_freq, sd_freq, phi, CH1, CP1, n_obs)
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
            CH1 = unique(CH1),
            CP1 = unique(CP1),
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
  #labs(title = expression(paste(c[H]^{(1)},"=",CH1_select," ",c[P]^{(1)},"=",CP1_select))) +
  labs(title = bquote(c[H]^{(1)}~"=" ~.(CH1_select) ~ c[P]^{(1)} ~"="~ .(CP1_select))) 

pdf(paste0(outprefix,"_sd_tint.pdf"), width = 13, height = 6.5)
print(sd_plt)
dev.off()

write_tsv(dpts, paste0(outprefix,"_interval_means_seed.tsv"))
write_tsv(sumover_investigation_dat,  paste0(outprefix,"_interval_means_overall.tsv"))
write_tsv(investigation_dat, paste0(outprefix, "_all_seeds_combined.tsv"))

pick_seed_10 <- sample(sum_investigation_dat |> pull(seed) |> unique(), 10, replace = FALSE)
pick_seed_2 <- sample(pick_seed_10, 2, replace = FALSE)

b <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

filter_investigation_dat <- filter_investigation_dat |>
  mutate(seed_level = factor(seed, levels=sort(unique(as.numeric(seed)))))


interval_img_plt <- ggplot(filter_investigation_dat |> filter(seed%in%pick_seed_10), aes(Genotype, interval_level)) +
  geom_tile(aes(fill = mean_freq)) +  
  geom_text(aes(label = round(mean_freq, 2)), color = "white", size = 2.2) +
  scale_fill_gradientn(limits = c(0,1),
                       colours=c("navyblue", "lightblue3", "lightblue", "darkorange1", "darkorange3", "brown3"),
                       breaks=b, labels=format(b)) +
  facet_grid(cols = vars(Species), rows = vars(seed), scales = "free_x") + 
  theme_bw() +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        title = element_text(size = 12, face = "bold")) +
  ylab("time") 

pdf(paste0(outprefix, "_interval_img.pdf"), width = 8, height = 10)
print(interval_img_plt)
dev.off()

# This plot is structured as follows:
# facet_wrap for genotypes with four columns and two rows
# each facet has the time intervals on the x-axis and the average frequency of a given genotype on the y-axis
# Color: The species of the presented genotype
# Datapoints for a given Genotype/Species/seed combination are connected with a geom_segment
geno_plt <- ggplot(dpts |> filter(interval_level != lastfirst), aes(interval_level, mean_freq)) +
  facet_wrap(vars(Genotype), ncol=4, nrow=2) +
  geom_path(aes(color=Species, group = interaction(Species, seed)), alpha=0.2) +
  #geom_jitter(aes(color=Species, shape = interval_level), height = 0, width = 0.2, alpha=0.7) + 
  geom_point(aes(color=Species), alpha=0.2) +
  scale_color_manual(values = c("#490092FF", "#006DDBFF"  ,"#004949FF")) +
  #geom_path(data = sumover_investigation_dat_test, aes(interval_level, mean_freq, group = Species, color = Species)) +
  #geom_point(data = sumover_investigation_dat_test, aes(interval_level, mean_freq, color=Species)) +
  #scale_shape_manual(values = c(16, 8, 17, 3, 15), name = "time interval") +
  #scale_alpha_discrete(range = c(0.2,1), guide = 'none') +
  theme_bw() + 
  xlab("Time interval") + 
  ylab("Mean genotype frequency") +
  guides(color = guide_legend(override.aes = list(alpha = 1)), alpha = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, size = 12)) +
  labs(title = bquote(c[H]^{(1)}~"=" ~.(CH1_select) ~ c[P]^{(1)} ~"="~ .(CP1_select))) 


pdf(paste0(outprefix,"_genotypes_seed.pdf"), width = 8.5, height = 5)
print(geno_plt)
dev.off()
  
#  Now summarize across seeds

outtest <- lapply(investigation_list, function(seed_dat){
  
  seed_dat_trans <- seed_dat |> 
    # Combine the Species and Genotype information into a single column speGeno
    mutate(speGeno = paste(Species, Genotype, sep="_"))
  
  #Pull all Species/Genotype combinations which got lost
  lost <- seed_dat_trans |> 
    filter(time == 100000 & freq == 0) |> 
    pull(speGeno)
  
  
  # Pull all Species/Genotype combinations which are still present
  present <-  seed_dat_trans |> 
    filter(time == 100000 & (freq > 0) ) |> 
    pull(speGeno)
  
  # Extract the last time point in the simulation 
  # where a genotype which was lost throughout the simulation is still present
  times_to_select <- sapply(lost, function(x){
    seed_dat_trans |> 
      filter(freq > 0 & speGeno == x) |>
      pull(time) |>
      max()
  })
  
  times_to_select_tib <- tibble(time=times_to_select, 
                                Genotype = gsub("[A-Z][0-2]*_","", names(times_to_select)),
                                Species = gsub("_[0-1]+$","",  names(times_to_select)),
                                seed = unique(seed_dat[["seed"]]),
                                speGeno = names(times_to_select))
  
  times_to_select_tib <- times_to_select_tib |>
    arrange(time) |>
    mutate(ext_order =1:n())
  
  dat_subset <- seed_dat_trans  |> 
    filter(time < (max(times_to_select)*1.5) & freq > 0)
  
  snaps_losses <- seed_dat_trans |>
    filter(time%in%c(0, as.numeric(times_to_select))) |>
    left_join(times_to_select_tib |> select(time, speGeno) |> rename("lost_geno"="speGeno"), by="time")
  
  seed_ext_heat_plt <- ggplot(snaps_losses, aes(Genotype, as.factor(time))) +
    geom_tile(aes(fill = freq)) +  
    geom_text(aes(label = round(freq, 6)), color = "white", size = 2.2) +
    scale_fill_gradientn(limits = c(0,1),
                         colours=c("navyblue", "lightblue3", "lightblue", "darkorange1", "darkorange3", "brown3"),
                         breaks=b, labels=format(b)) +
    facet_wrap(vars(Species), ncol =3, scales = "free_x") + 
    theme_minimal() +
    theme(strip.text.x = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          legend.title = element_blank(),
          title = element_text(size = 12, face = "bold")) +
    ylab("time") 

  return(list(times_to_select_tib, dat_subset, snaps_losses, seed_ext_heat_plt))
  
})

datapts <- lapply(outtest, function(x){
  return(x[[2]])
})

extpts <- lapply(outtest, function(x){
  return(x[[1]])
})

losspts <- lapply(outtest, function(x){
  return(x[[3]])
})

ls_image_list <- lapply(outtest, function(x){
  return(x[[4]])
})

pdf(paste0(outprefix,"_seed_ext_freqs_heat.pdf"), width = 12 , height =5)
for(i in 1:length(ls_image_list)){
  print(ls_image_list[[i]])
}
dev.off()

datpts_dat <- bind_rows(datapts)


extpts_dat <- bind_rows(extpts) |>
  mutate(Genotype = factor(Genotype, levels = c("000","100","010","001","110","101","011","111")))

losspts_dat <- bind_rows(losspts) |>
  mutate(Genotype = factor(Genotype, levels = c("000","100","010","001","110","101","011","111")))



extdyn1_plt <- ggplot(datpts_dat |> filter(seed%in%pick_seed_10[1:5]), aes(time, freq)) +
  geom_line(aes(color=Genotype)) +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  geom_point(data= (datpts_dat |> filter(time == 0 & seed%in%pick_seed_10[1:5])), aes(time, freq, color=Genotype), size=3, alpha = 0.6) +
  facet_grid(rows = vars(Species), cols = vars(seed), scales = "free_x") +
  theme_bw() +
  geom_vline(data=extpts_dat |> filter(seed%in%pick_seed_10[1:5]), aes(xintercept=time,color=Genotype, linetype = Genotype), 
             key_glyph = "path", linewidth = 1)  +
  scale_shape_manual(values=c(0,2,6,1,3,7,8,9), breaks = c("000","100","010","001","110","101","011","111"))

extdyn2_plt <- ggplot(datpts_dat |> filter(seed%in%pick_seed_10[6:10]), aes(time, freq)) +
  geom_path(aes(color=Genotype)) +
  #geom_point(aes(color=Genotype, shape = Genotype), alpha = 0.8) +
  scale_color_manual(values =  c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  geom_point(data= (datpts_dat |> filter(time == 0 & seed%in%pick_seed_10[6:10])), aes(time, freq, color=Genotype), size=3, alpha = 0.6) +
  facet_grid(rows = vars(Species), cols = vars(seed), scales = "free_x") +
  theme_bw() +
  geom_vline(data=extpts_dat |> filter(seed%in%pick_seed_10[6:10]), aes(xintercept=time,color=Genotype, linetype = Genotype), 
             key_glyph = "path", linewidth = 1) +
  scale_shape_manual(values=c(0,2,6,1,3,7,8,9), breaks = c("000","100","010","001","110","101","011","111"))

pdf(paste0(outprefix,"_extdyn.pdf"), width = 11, height = 5)
print(extdyn1_plt)
print(extdyn2_plt)
dev.off()

ext_sum <- extpts_dat |> 
  group_by(Species, Genotype) |>
  summarize(n=n())

extbxplt <- ggplot(extpts_dat, aes(time, Genotype)) +
  facet_grid(rows=vars(Species)) + 
  geom_boxplot(aes(color=Genotype), alpha = 0.8) +
  geom_point(aes(color=Genotype)) +
  scale_color_manual(values = c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  theme_bw() + 
  geom_text(data=ext_sum, aes(Inf, Genotype, label = paste("n=",n)), hjust = 1.2)



pdf(paste0(outprefix,"_extbxplt.pdf"), width = 7, height = 5)
print(extbxplt)
dev.off()

extorder_seed_plt <- ggplot(extpts_dat |> arrange(time), aes(time, as.factor(seed))) +
  #geom_path(aes(group=speGeno, linetype = speGeno)) +
  geom_point(aes(color = Species, shape = Genotype), size = 3) +
  scale_color_manual(values = c("#490092FF", "#006DDBFF"  ,"#004949FF")) +
  scale_shape_manual(values = c(0:3,5:8)) + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(linewidth = 2),
        panel.grid.major.x = element_blank()) +
  ylab("Seed")


pdf(paste0(outprefix,"_extbxplt_seed.pdf"), width = 7.5, height = 4.5)
print(extorder_seed_plt)
dev.off()

extorder_plt <- ggplot(extpts_dat, aes(ext_order, as.factor(seed))) +
  geom_point(aes(shape = Species, color = Genotype), size = 3) +
  scale_color_manual(values = c("000"="gray80", "100"="#D9565CFF","010"="#F28A8AFF",
                                 "001"="brown4","110"="#1BB6AFFF","101"="#088BBEFF","011"="#172869FF","111"="gray20")) +
  theme_bw() +
  scale_shape_manual(values = c(16, 17, 15)) +
  theme(panel.grid.major.y = element_line(linewidth = 2),
        panel.grid.major.x = element_blank()) +
  ylab("Seed") +
  xlab("Extinction order") +
  scale_x_continuous(breaks = 1:max(extpts_dat[["ext_order"]])) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 12, vjust = -1))

pdf(paste0(outprefix, "_extinction_order.pdf"), width = 4.5, height = 8)
print(extorder_plt)
dev.off()

write_tsv(extpts_dat, paste0(outprefix, "_extpts_dat.tsv"))
write_tsv(losspts_dat, paste0(outprefix, "_losspts_dat.tsv"))


pplt_final <- investigation_dat |> 
  filter(Species == "P" & time >  95000 & freq > 0 & seed%in%pick_seed_10)

pplt_check <- pplt_final |> 
  group_by(seed, Species, Genotype) |> 
  summarize(n=n()) |>
  pull(n)

if(any(pplt_check < 10)){
  pplt_final <- investigation_dat |> 
    filter(Species == "P" & seed%in%pick_seed_10) |>
    group_by(seed, Species, Genotype) |>
    slice_tail(n=5)
}

pat_dynamics_plt <- ggplot(pplt_final |> filter(freq > 0), aes(time, freq)) +
  geom_line(aes(color=as.factor(seed)))  +
  theme_bw() +
  scale_color_paletteer_d("ggsci::nrc_npg", name = "seed") +
  facet_wrap(vars(Genotype), nrow = 2) + 
  ylim(0,1)

pdf(paste0(outprefix, "_pat_end.pdf"), width = 8.5, height = 6)
print(pat_dynamics_plt)
dev.off()


hostplt_final <- investigation_dat |> 
  filter(Species %in% c("H1","H2") & time >  95000 & freq > 0 & seed%in%pick_seed_10)


if(any(pplt_check < 10)){
  hostplt_final <- investigation_dat |> 
    filter(Species %in% c("H1","H2") & seed%in%pick_seed_10) |>
    group_by(seed, Species, Genotype) |>
    slice_tail(n=5)
}


host_dynamics_plt <- ggplot(hostplt_final |> filter(freq > 0), aes(time, freq)) +
  geom_line(aes(color=as.factor(seed)))  +
  theme_bw() +
  scale_color_paletteer_d("ggsci::nrc_npg", name = "seed") +
  facet_grid(rows=vars(broader_cat), cols = vars(Species))+ 
  ylim(0,1)

pdf(paste0(outprefix, "_host_end.pdf"), width = 8.5, height = 10)
print(host_dynamics_plt)
dev.off()

# Calculate the mean for the last 5,000 time steps
if(!any(pplt_check < 10)){
  mean_data <- investigation_dat |> 
    filter(time > 95000 ) %>%
    group_by(seed, Species, Genotype) %>%
    summarize(mean_freq = mean(freq, na.rm = TRUE), .groups = "drop")
}else{
  mean_data <- investigation_dat |> 
    group_by(seed, Species, Genotype) |>
    slice_tail(n = 5) |>
    summarize(mean_freq = mean(freq, na.rm = TRUE), .groups = "drop") 
}

mean_data_last_wider <- mean_data |> 
  mutate(mean_freq = round(mean_freq, digits=2)) |> 
  pivot_wider(values_from = mean_freq, names_from = seed, names_prefix = "seed_")

if(!any(pplt_check < 10)){
  mean_data_last_wider <- mean_data_last_wider |>
    mutate(cat = "last_5000")
  pplt_final <- pplt_final |>
    mutate(cat = "last_5000")
  hostplt_final <- hostplt_final |>
    mutate(cat = "last_5000")
}else{
  mean_data_last_wider <- mean_data_last_wider |>
    mutate(cat = "last_5")
  pplt_final <- pplt_final |>
    mutate(cat = "last_5")
  hostplt_final <- hostplt_final |>
    mutate(cat = "last_5")
}

write_tsv(mean_data_last_wider, paste0(outprefix,"_last5000_means.tsv"))
write_tsv(hostplt_final, paste0(outprefix,"_last_pointsH.tsv"))
write_tsv(pplt_final, paste0(outprefix,"_last_pointsP.tsv"))

mean_data <- mean_data |>
  mutate(broader_cat = case_when(
    Species%in%c("H1","H2") & Genotype%in%c("100","010") ~ "R_priv(100/010)",
    Species%in%c("H1","H2") & Genotype%in%c("101","011") ~ "R_double(101/011)",
    TRUE ~ Genotype
  )) |>
  mutate(broader_cat = factor(broader_cat, levels = c("000","100","010","001","R_priv(100/010)","101","011","110","R_double(101/011)","111")))

if(!any(pplt_check < 10)){
  average_dyn_hend_plt <- ggplot(investigation_dat |> filter(time > 95000 & Species%in%c("H1","H2") & freq > 0 & seed%in%pick_seed_10),
         aes(time, freq)) +
    geom_line(data = (investigation_dat |> filter(time > 95000 & Species%in%c("H1","H2") & freq > 0 & !seed%in%pick_seed_10)), linewidth = 0.15, color = "gray")  +
    geom_line(data =investigation_dat |> filter(time > 95000 & Species%in%c("H1","H2") & freq > 0 & seed%in%pick_seed_10), aes(color= as.factor(seed)), linewidth = 0.2)  +
    geom_hline(data = (mean_data |> filter(Species%in%c("H1","H2") & mean_freq > 0 & seed%in%pick_seed_10)), aes(yintercept = mean_freq), 
               linetype = "solid", linewidth = 0.8, alpha = 0.7, color =  "gray") +
    geom_hline(data = (mean_data |> filter(Species%in%c("H1","H2") & mean_freq > 0 & seed%in%pick_seed_10)), aes(yintercept = mean_freq, color = as.factor(seed)), 
                linetype = "solid", linewidth = 1, alpha = 0.7) +
    theme_bw() +
    scale_color_paletteer_d("ggsci::nrc_npg", name = "seed") +
    facet_grid(cols = vars(Species), rows = vars(broader_cat)) + 
    ylim(0,1)
  
  
  pdf(paste0(outprefix,"_average_dyn_hend.pdf"), height = 5, width = 8)
  print(average_dyn_hend_plt)
  dev.off()
  
  
  average_dyn_pend_plt <- ggplot(investigation_dat |> filter(time > 95000 & Species%in%c("P") & freq > 0 & seed%in%pick_seed_10),
                                 aes(time, freq)) +
    geom_line(data = (investigation_dat |> filter(time > 95000 & Species%in%c("P") & freq > 0 & !seed%in%pick_seed_10)), linewidth = 0.15, color = "gray")  +
    geom_line(data =investigation_dat |> filter(time > 95000 & Species%in%c("P") & freq > 0 & seed%in%pick_seed_10), aes(color= as.factor(seed)), linewidth = 0.2)  +
    geom_hline(data = (mean_data |> filter(Species%in%c("P") & mean_freq > 0 & seed%in%pick_seed_10)), aes(yintercept = mean_freq), 
               linetype = "solid", linewidth = 0.8, alpha = 0.7, color =  "gray") +
    geom_hline(data = (mean_data |> filter(Species%in%c("P") & mean_freq > 0 & seed%in%pick_seed_10)), aes(yintercept = mean_freq, color = as.factor(seed)), 
               linetype = "solid", linewidth = 1, alpha = 0.7) +
    theme_bw() +
    scale_color_paletteer_d("ggsci::nrc_npg", name = "seed") +
    facet_wrap(vars(Genotype), ncol = 4, nrow = 2) +
    ylim(0,1)
  
  
  pdf(paste0(outprefix,"_average_dyn_pend.pdf"), height = 5, width = 10)
  print(average_dyn_pend_plt)
  dev.off()
}

initial_cond <- investigation_dat |> 
  filter(time == 0) |>
  arrange(seed) |>
  select(seed, Species, Genotype, freq) 


initial_plt <- ggplot(initial_cond, aes(Genotype, as.factor(seed))) +
  geom_tile(aes(fill = freq)) +  
  geom_text(aes(label = round(freq, 2))) +
  scale_fill_gradientn(limits = c(0,1),
                       colours=c("navyblue", "lightblue3", "lightblue", "darkorange1", "darkorange3", "brown3"),
                       breaks=b, labels=format(b)) +
  facet_wrap(vars(Species), ncol =3, scales = "free_x") + 
  theme_minimal() +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        title = element_text(size = 12, face = "bold")) +
  ylab("Seed") +
  ggtitle("Initial frequencies")
  


pdf(paste0(outprefix,"_initial.pdf"), width = 10, height = 4.5)
print(initial_plt)
dev.off()
  
initial_cond_wider <- initial_cond |>
  mutate(freq = round(freq, digits = 2)) |> 
  pivot_wider(names_from = seed, values_from = freq) 

write_tsv(initial_cond_wider, paste0(outprefix,"_initial_cond.tsv"))

