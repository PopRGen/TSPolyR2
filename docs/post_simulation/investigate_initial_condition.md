i# investigate_initial_condition.R

## Overview 

This script can be used to produce several diagnostic plots summarizing differences for several simulations which have been run with the exact same parameters, but different initial conditions

## Dependencies

- tidyverse
- argparser
- paletteer
- ggplot2
- patchwork
- gridExtra

## Usage

```{bash,eval=F}
Rscript --vanilla investigate_initial_condition.R [OPTIONS]
```

## Options

- `--prefix`: Prefix to append to all simulations. Default: testrun_small. Type: character.
- `--suffix`: Suffix to append to all filenames. Default="2024-11-15. Type: character. 
- `--cH_2`: Value of cH_2 to use for the simulations. Default: 3. Type: numeric.
- `--cP_2`: Value of cP_2 to use for simulations. Default: 0. Type: numeric.
- `--indir`: Top level directory for the run. Default: "/scratch/hm2840/modelling_rerun/fsched_small_rand_cleanup_ext_rerun/testrun_small_CH2_3_CP2_0_2024-11-15/all_dynamics_long. Type: character.
- `--plotdir`: Directory for storing the simulations. Default: Figures_seed_comparison/extended_analyses. Type: character.
- `--tselect`: Time from which onwards frequencies should be averaged. Default:  40000. Type: numeric.
- `--phi`: Value of phi to analyze. Default:  0.5. Type: numeric
- `--cH_1`: Value of cH_1 to analyze. Default:  0.1. Type: numeric
- `--cP_1`: Value of cP_1 to analyze. Default:  0.2. Type: numeric
- `--nrep`: number of different initial seed. Default: 50. Type: numeric.
- `--seed_start`: lower seed to start with. Default: 400. Type: numeric
- `--seed_interval`: seed interval tp ise. Default: 400. Type: numeric.

## Output

- `<outprefix>_sd_tint.pdf`:  Average frequency of each genotype in a given time interval. Facet plot. Columns are time intervals, rows are species. Genotypes are on the x-axis. The result for each simulation is plotted as point. y-axis average frequency of the genotype in the inveral.

- `<outprefix>_interval_means_seed.tsv`: File contains the mean genotype frequency and standard deviation of a given genotype for the time intervals 0-20,000; 20,000 - 40,000; 40,000 - 80,000; 80,000 - 100,000, 100,000 for each seed.

    | Column name | Description |
    | :--- | :--- |
    | interval_level | For which time interval all time points are assigned to 20,000 t intervals |	
    | seed | seed for the initial genotype frequencies |
    | Species | species |
    | Genotype | genotype in the species |
    | mean_freq | mean frequency of the genotype in the time interval |
    | sd_freq | standard deviation of the genotype frequency in the time interval |
    | phi | value of $\phi$ for the simulation |
    | CH1 | value of $c_H^{(1)} for the simulation |
    | CP1 | value of $c_P^{(1)}$ for the analyzed simulation |
    | n_obs | number of observations. I think currently not calculated correctly |

- `<outprefix_interval_means_overall.tsv`: File contains the mean genotype frequency and standard deviation of a given genotype for the time intervals 0-20,000; 20,000 - 40,000; 40,000 - 80,000; 80,000 - 100,000, 100,000 across all seeds.

    | Column name | Description |
    | :--- | :--- |
    | interval_level | time interval analyzed |	
    | Species | species |
    | Genotype | genotype |
    | mean_freq | mean frequency of the genotype across seeds |
    | sd_freq | standard deviation of the genotype frequency across all seeds |
    | phi | value of $\phi$ for the simulation |
    | CH1 | cost of resistance for the simulations ($c_H^{(1)}$) |
    | CP1 | cost of virulence for the simulations ($c_P^{(1)}$) |	
    | n_obs | number of observations. I think currently not calculated correctly | 

- `<outprefix>_all_seeds_combined.tsv`: Contains the combined dynamics across all seeds in long format.

    | Column | Description|
    | :--- | :--- |	
    | Time | time point|
    | Type | Short form form for the genotype information. I.e. H1_0 |	
    | Species | Species |
    | ID | Short ID of the genotype.|
    | freq | frequency of the genotype |	
    | Genotype | decoded genotype (allelic states) |	
    | cumsum_freq | cumulative sum of all genotype frequency in the species including the genotype currently analyzed |
    | lower_freq | cumulative sum of all genotype frequency in the species before the current genotype is analyzed |
    | phi | value of $\phi$ used for the simulation|	
    | CH1 | cost of resistance ($c_H^{(1)}$) |	
    | CP1 | cost of virulence ($c_P^{(1)}$) |	
    | CH2 | shape of the host cost function ($c_H^{(2)}$)|
    | CP2 | shape of the pathogen cost function ($c_P^{(2)}$)|
    | seed | seed for intializing the genotype frequencies for the given simulation|
    | interval | coarse grained time interval this data point is part of |
    | interval2 | fine-grained time interval this data point is part of|
    | interval_level | see interval (this was the column formatted as factor)|
    | interval2_level |see interval2 (this was the column formatted as factor)|	
    | broader_cat | broader category, mainly used to describe similar genotypes between both host species such as R_priv(100/010) |

- `<outprefix>_interval_img.pdf`: Image plot of the average frequencies of all genotypes in 20,000 time span interval. Facet grid plot with three columns (species) and 10 rows (10 different seeds). Each facet corresponds to an image plot of the average frequency of a given gentoype (x-axis) in a given time interval.

- `<outprefix>_genotypes_seed.pdf`: Average of all genotypes for different 20,000 time span intervals. Facet wrap plot for all potential genotypes (2 rows x 4 columns). First row (000, 100, 010, 001), second row (110,101,011,111). Each facet show the average frequency of the given genotype (y-axis) in a given time interval. Results are shown for all seeds and and species and colored by species. Results for a single seeds for a single species are connected by a line. Each of this lines has some alpha value. The average across seeds is drawn with alpha = 1. 

- `<outprefix>_extdyn.pdf`: Plot of the extinction dynamics for ten seeds. Shown are to facet grid plots. Each of them has five columns (five different seeds) and three rows (three species). Time points were a given genotype got extinct in a given species are indicated by a vertical bar.  

- `<outprefix>_extbxplt.pdf`: Boxplot of the extinction order for all different initial genotype frequencies. Facet plot with three rows (bottom: host H, middle: host M, bottom: pathogen). Shown are the exinction times (x-axis) of single genotype (y-axis) in the species and the corresponding boxplot for all $nreps$ simulations.

- `<outprefix>_extbxplt_seed.pdf`: Extinction order for each seed. Extinction times are on the x-axis and the seeds on the y-axis. Different species are indicated by color and the genotypes by shape.

- `<outprefix>_extinction_order.pdf`: Contains the order of extinctions for all seeds. Extinction ranks are on the x-axis, the different seeds on the y-axis. The extinct genotype is indicated by the color and its species by the shape (circle host $H$, triangle host $M$, square pathogen).

- `<outprefix>_extpts_dat.tsv`: Extinction time of the genotypes which got extinct across all seeds.

    | Column | Description |
    | :--- | :--- |
    | time | extinction time |
    | Genotype | genotype that got extinct |	
    | Species | species the genotype belongs to |
    | seed | seed of the simulation |
    | speGeno | combination of the species and the genotype |
    | ext_order | the how maniest genotype was this genotype to get extinct in the simulation?

- `<outprefix>_losspts_dat.tsv`: Frequencies of all genotypes at the time point before a genotype gets extinct. Information from all seeds has been put together into a single file.

    | Column | Description |
    | :--- | :--- |
    | time | time point immediately before one or several genotypes got extinct |
    | Type | the genotype in the H1_000 ... notation |
    | Species | the species of the genotype |
    | ID | the numeric ID of the genotype in the species |
    | freq | frequency of the genotype |
    | Genotype | the genotype in  000 notation | 	
    | cumsum_freq | cumulative frequency including this genotype in the given simulation at the given time point |
    | lower_freq | the cumulative frequency of all genotypes up to the current genotype in the given simulation at the given time point |
    | phi | the value of $\phi$ for the given simulation |
    | CH1 | cost of resistance ($c_H^{(1)}$) |	
    | CP1 | cost of virulence ($c_P^{(1)}$) |	
    | CH2 | shape of the host cost function ($c_H^{(2)}$)|
    | CP2 | shape of the pathogen cost function ($c_P^{(2)}$)|
    | seed | the seed for the given simulation |
    | speGeno | the species genotype combination |
    | lost_geno | the genotype being lost in the nexrt time step.

- `<outprefix>_pat_end.pdf`: Frequency dynamics (y-axis) of all maintained pathogen genotypes for all data points in the time interval 95,000 - 100,000 (x-axis) across all seeds. Each seed is indicated by a different color.

- `<outprefix>_host_end.pdf`: Frequency dynamics of all maintained host H genotypes (left) amd host M genotypes (right). There is one panel for each genotypes. The dynamics for different seeds are indicated by color.

- `<outprefix>_last5000_means.tsv`: Mean frequency of all maintained genotypes for interval 95,000 - 100,000 across all seeds for the given combination. One row for each gentoype. One column for each seed. Indicated are the mean frequencies across all data points in the time interval 95,000 - 100,000.

- `<outprefix>_last_pointsH.tsv`: Potentially no longer required.

- `<outprefix>_last_pointsP.tsv`: Potentiall no longer required.

- `<outprefix>_average_dyn_hend.pdf`: Facet plot of the frequency dynamics (y-axis) for 10 different initial conditions (color) during time 95,000 - 100,000 (x-axis) for all maintained genotypes (rows) in the two hosts (columns).

- `<outprefix>_average_dyn_pend.pdf`: Facet plot of the frequency dynamics (y-axis) for 10 different initial conditions (color ) during time 95,000 - 100,000 (x-axis) for all maintained genotypes (rows). Note for consistency the same simulations are chosen as in previous plot.

- `<outprefix>_initial.pdf`: Image plots of the initial frequencies across seeds. Three image plots. Left host $H$, middle host $M$ and right pathogen. Each image plot has the genotypes on the x-axis and the seeds on the y-axis.

- `<outprefix>_initial_cond.tsv`: Overview of all initial conditions for the given parameter combination. One row per genotype. One column per seed. Values correspond to the rounded intial frequencies.

