---
title: "README Two host/generalist pathogen model"
author: "Hanna Maerkle"
date: "12/10/2025"
output: html_document
---


# TSPolyR2


This is the updated README to simulate interactions between to closely related hosts and a generalist pathogen. The model is based on the non-ecological model in Ashby and Boots (2017) [doi: 10.1111/ele.12734](10.1111/ele.12734).
In short we extended the model to incorporate two phylogenetically related host species interacting with a generalist pathogen. 

The main parameters which have been investigated in the model are:



## Cloning the repository and compiling the executable:

1. Clone the git repository.  

    ```{bash, eval=FALSE}
    git clone https://github.com/PopRGen/TSPolyR2.git
    ```

2. Run the make file to compile the exectuable.  

    ```{bash, eval=FALSE}
    cd TSPolyR2/src
    make
    ```

## Running a single simulation

    ```{bash}
    ./Non_ecological_clean
    ```

## Running simulations for a fixed combination of shape parameters ($\Xi_H$=$\Xi_M$, $\Xi_P$) with random values for the maximum costs ($\Omega_H$=$\Omega_M$, $\Omega_P$) and random inital genotype frequencies

For instructions see [docs/running_simulations/simulations_random.R](./docs/running_simulations/simulations_random.R).


### Summarize the results

The results can be summarized with [post_simulation/summarize_analyze_random.R](./docs/post_simulation/summarize_analyze_random.R).

A slurm-submission script to summarize the simulations for $\Xi_H=\Xi_M=\Xi_P=3$, $\phi=0.5$ and RSEED=1700 can be found [here](./docs/post_simulation/run_summary_random.md).

### Producing the main figures in the manuscript

All files for producing the figures in the manuscript can be found in the directory [post_simulation](post_simulation).



## Running simulations for fixed combinations of $c_P^{(1)}$, $c_H^{(1)}$, $c_H^{(2)}$, $c_P^{(2)}$ and $\phi$ with random intial frequencies.

Follow steps 1. and 2. above. Use the script [src/Run_simulation.R] to run a simulation inside an array job. 
Specific instructions and the corresponding slurm submission script can be found here: [docs/running_simulations/simulations_repeated.md](./docs/running_simulations/simulations_repeated.md)
 
















