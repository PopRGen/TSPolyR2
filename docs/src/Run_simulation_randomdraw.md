<style type="text/css">
  body{
  font-size: 8pt;
}
</style>

# Run_simulation_randomdraw.md

## Overview

Script to run the model `nreps` times by randomly drawing values for $c_H^{(1)}$, $c_P^{(1)}$ and $\phi$ for each simulation, while keeping all other parameters fixed. Note: If the upper and lower limit are defined to be same value, the given parameter will be fixed among all `nreps` simulations. The script automatically draws a random seed for each of the `nreps` simulations to intialize the random number generator in the C++-program which is used to simulate the dynamics.

## Dependencies

* tidyverse

* argparser

## Usage

```
Rscript --vanilla Run_simulation_randomdraw.R \
    --cH_2 3 \
    --cP_2 -3 \
    --executable TSPolyR/src/Non_ecological_clean  \
    --cH1_min 0.01 \
    --cH1_max 0.3 \
    --cP1_min 0.06 \
    --cP1_max 0.2 \
    --phi_min 0.5 \
    --phi_max 0.5 \
    --init r \
    --excheck 1 \
    --nreps 100 
```

This would run 100 simulations (`--nreps 100`) for $c_H^{(2)}=3$ (`--cH_2 3`)  and $c_P^{(2)}=-3$ (`--cP_2 -3`) and $\phi=0.5$ (as `phi_min 0.5`, `--phi_max 0.5`). For each simulation genotype frequencies are initialized at random (`--init r`). A value for $c_H^{(1)}$ is drawn from a uniform distribution $\mathcal(U)(0.01,0.3)$ (`--cH1_min 0.01`, `--cH1_max 0.3`) and a value for $c_P^{(1)}$ is drawn from a uniform distribution $\mathcal(U)(0.06,0.2)$ (`--cP1_min 0.06`, `--cP1_max 0.2`). Random initialization of genotype frequencies is performed (`--init r`). During the simulations recurrent checks are performed if genotypes got extinct (`--excheck 1`). 


## Options

- `--prefix`: Prefix to append to all simulations. Default: `Non_ecological`. Type: character.

- `--suffix`: Suffix to append to all filenames. Default: `format(Sys.Date(),"%Y-%m-%d")`. Type: character.

- `--executable`: Location of the compiled C++-code used for running the simulations. Default: `Non_ecological_clean`. Type: character.

- `--cH_2`:  Value of cH_2 to use for the simulations. Default: 3. Type: character.

- `--cP_2`: Value of cP_2 to use for the simulations. Default: 3. Type: character.

- `--block1`: Which locus to block of in host1. Note is zero-based. Default: 0. Type: numeric.

- `--block2`: Which locus to block of in host2. Note is zero-based. Default: 1. Type: numeric.

- `--sigma`: Value of $\sigma$ to use in the simulations. Default: 0.85. Type: numeric.

- `--betaH`: Value of $\beta_H$ to use in the simulations. Default: 1. Type: numeric.

- `--betaP`: Value of $\beta_H$ to use in the simulations. Default: 1. Type: numeric.

- `--nloci`: Number of loci to simulate. Default:3. Type: numeric.

- `--init`: How to initialize the genotype frequencies. Default: 'f'. Type: character.

- `--excheck`: Check if genotypes went extinct. 0= no, 1 = yes. Default: 0. Type: numeric.

- `--cH1_min`: Lower limit of the uniform distribution from which the $c_H^{(1)}$ are drawn. Default: 0.1. Type: numeric.

- `--cH1_max`: Upper limit of the uniform distribution from which the $c_H^{(1)}$ are drawn. Default: 0.3. Type: numeric.

- `--cP1_min`: Lower limit of the uniform distribution from which the $c_P^{(1)}$ are drawn. Default: 0.1. Type: numeric.

- `--phi_min`: Lower limit of the uniform distribution from which the $\phi$ are drawn. Default: 0.5. Type: numeric.

- `--phi_max`: Lower limit of the uniform distribution from which the $\phi$ are drawn. Note: if `phi_min`=`phi_max`, simulations are for a fixed value of phi. Default: 0.5. Type: numeric.

- `--rseed`: seed to use the RNG in R. Default: 2. Type: numeric.

- `--outdir`: Path to the output directory in which the results should be stored. Default: `../../results_random`. Type: character.

- `--nreps`: Number of simulations to run. Default: 10. Type: numeric.


## Additional information

## Output

* `<outdir>/<prefix>_CH2_<cH_2>_CP2_<cP_2>_random_combis_<suffix>/simulations`. A directory with the results for all `nreps` simulations.

    * `CH2_<cH_2>_CP2_<cP_2>_random_r_<nreps>_RSEED_<rseed>.tsv`: tab-separate file with four columns. Contains the values of $c_H^{(1)}$ (column `cH_1`), $c_P^{(1)}$ (column `cP_1`), $\phi$ (column `phi`) and the seed (column `seed`) to intialize the random number generator of the C++ program for each of the `nreps` simulations

    * `CH2_<cH_2>_CP2_<cP_2>_random_r_<nreps>_RSEED_<rseed>.tar`: A tar archive which contains the result for each single simulation. Archive is automatically created to not hit file number limitations.

    * `CH2_<cH_2>_CP2_<cP_2>_random_r_<nreps>_RSEED_<rseed>_end.txt`: File with the frequencies at the end of each simulation. 
