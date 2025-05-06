
# Plot_fluctuations.R

## Overview

This script to plot the dynamics for a single simulation.

## Dependencies

- tidyverse
- argparser
- paletteer
- ggplot2
- patchwork
- ggpubr

## Usage

```{bash}
Rscript --vanilla Plot_fluctuations.R [OPTIONS]
```

## Options

- `--prefix`: Prefix which have been appended to all simulations of the set. Default: ext_reps. Type: character.
-  `--suffix`: Suffix to append to all filenames. Default: 2025-01-07. Type: character.
- `--phi`: Value of phi to use for the simulations. Default: 0.5, Type: numeric.
- `--cH_2`: Value of cH_2 to use for the simulations. Default: 3. Type: numeric.
- `--cP_2`: Value of cP_2 to use for simulations. Default:3. Type: numeric.
- `--cH_1`: Value of cH_1 to use for the simulations. Default: 0.1. Type: numeric.
- `--cP_1`: Value of cP_1 to use for simulations. Default: 0.1. Type: numeric.
- `--parent_outdir`: Top level directory for results. Default: results_r50. Type: character.
- `--seed`: Seed to use for the simulation. Default: 1600. Type: numeric.
- `--plotdir`: Directory for storing the simulations. Default: Figures. Type: character.


## Output

- `<plotdir>/dynamics_cH2_<CH2>_cP2_<CP2>_phi_<phi>_cH1_<CH1>_cP1_<CP1>.pdf`: Plot of the dynamics over time.

- `Fig1.pdf`: Figure 1 as shown in the current version of the paper (assuming standard settings are used).
