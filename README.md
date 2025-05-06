---
title: "README Two host/generalist pathogen model"
author: "Hanna Maerkle"
date: "3/31/2025"
output: html_document
---


# TSPolyR2


This is the README to simulate interactions between to closely related hosts and a generalist pathogen. The model is based on the non-ecological model from Ashby and Boots (2017).


## Recipe 1: Single simulation:

1. Clone the git repository.  

    ```{bash, eval=FALSE}
    git clone https://github.com/PopRGen/TSPolyR.git
    ```

2. Run the make file to compile the exectuable.  

    ```{bash, eval=FALSE}
    cd TSPolyR/src
    make
    ```

3. Running a single simulation

    ```{bash}
    ./Non_ecological_clean
    ```


## Recipe 2: Running simulations for random combinations of $c_P^{(1)}$, $c_H^{(1)}$ for fixed values of $c_H^{(2)}$, $c_P^{(2)}$ and $\phi$.

Follow steps 1. and 2. from recipe 1.

Use the script [src/Run_simulation_randomdraw.R](./docs/src/Run_simulation_randomdraw.md) to run a fixed number of simulations drawing random values for $c_H^{(1)}$, $c_P^{(1)}$ and random initial genotype frequencies.

```{bash,eval=F}
cd src
Rscript --vanilla  TSPolyR2/src/Run_simulation_randomdraw.R \
    --cH_2 3 \
    --cP_2 3 \
    --cH1_min 0.05 \
    --cH1_max 0.3 \
    --cP1_min 0.01 \
    --cP1_max 0.3 \
    --phi_min 0.5 \
    --phi_max 0.5 \
    --init r \
    --rseed 15 \
    --nreps 100
```
This runs 100 simulations. For each simulation $c_H^{(2)} = c_P^{(2)} = 3$ and $\phi=0.5$. The values of $c_H^{(1)}$ for each simulation is drawn from $\mathcal{U}(0.01,0.3)$ and the value of $c_P^{(1)}$ for each simulation is drawn from $\mathcal{U}(0.01,0.3)$.

## Running simulations for fixed combinations of $c_P^{(1)}$, $c_H^{(1)}$, $c_H^{(2)}$, $c_P^{(2)}$ and $\phi$ with random intial frequencies.

Follow steps 1. and 2. from recipe 1. Use the script [src/Run_simulation.R](./docs/src/Run_simulation.md) to run a simulation.

```{bash,eval=F}
cd src 
Rscript --vanilla Run_simulations.R \
    --cH_2 3 \
    --cP_2 3 \
    --executable `pwd`/TSPolyR2/src/Non_ecological_clean \
    --phi 0.7 \
    --suffix 2025-01-27 \
    --parent_outdir results_r50 \
    --block1 0 \
    --block2 1 \
    --init r \
    --sigma 0.85 \
    --betaH 1 \
    --betaP 1 \
    --prefix ext_reps \
    --seed 16800 \
    --excheck 1 \
    --gsize small
```






















