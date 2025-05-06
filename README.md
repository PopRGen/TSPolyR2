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

Follow steps 1. and 2. from recipe 1. Use the script [src/Run_simulation.R](./docs/src/Run_simulation.md) to run a simulation inside an array job. Below is the script to run all pairwise combinations of $c_H^{(1)} \in \{0.1, 0.2, 0.3\}$, $c_P^{(1)} \in \{0.1, 0.2, 0.3\}$, for $c_H^{(2)}=c_P^{(2)}=3$ and $\phi=0.5$.  A copy of the corresponding simulation schedule `schedules_r50/schedule_ext_reps_CH2_3_CP2_3_phi_0.5_2025-01-07_random.txt` can be found in `auxillaries/schedules_r50`. Make sure that the submission script from below is placed into a parent directory, which contains a copy of this respository as a subdirectory and a `schedules_r50` subdirectory with the simulation schedule. Otherwise the submission script needs to be modified. 

```{bash,eval=F}
#!/bin/bash

#SBATCH --mem=2GB
#SBATCH --cpus-per-task=2
#SBATCH --time=10:00
#SBATCH --job-name=r50_CH2_3_CP2_3_phi_0.5
#SBATCH --output=logs/r50_rand50_CH2_3_CP2_3_phi_0.5_%A_%a.out
#SBATCH --error=logs/r50_rand50_CH2_3_CP2_3_phi_0.5_%A_%a.err
#SBATCH --tasks=1
#SBATCH --array=1-50

module load gcc/10.2.0
module load r/gcc/4.4.0

# Set the desired parameters for the simulation
SIMSCHEDULE="schedules_r50/schedule_ext_reps_CH2_3_CP2_3_phi_0.5_2025-01-07_random.txt"
BETAP=1
BETAH=1
SIGMA=0.85
BLOCK1=0
BLOCK2=1
EXCHECK=1

# These two paths are relative to the working_dir
WORKINGDIR=`pwd`
RESULTSDIR="${WORKINGDIR}/results_r50"
SIMULATIONDIR="${WORKINGDIR}/TSPolyR2/src"

ID=${SLURM_ARRAY_TASK_ID}

echo "WORKINGDIR:   ${WORKINGDIR}"
echo "SIMSCHEDULE:  ${SIMSCHEDULE}"
echo "BETAP:    ${BETAP}" 
echo "BETAH:    ${BETAH}" 
echo "SIGMA:    ${SIGMA}"
echo "BLOCK1:   ${BLOCK1}"
echo "BLOCK2:   ${BLOCK2}"
echo "RESULTSDIR:   ${RESULTSDIR}"
echo "SIMULATIONDIR:    ${SIMULATIONDIR}"



# Read params for current bunch of simulations
PARAMS=`sed -n "${ID} p" ${SIMSCHEDULE}`
parArr=(${PARAMS}) 

EXECUTABLE=${parArr[0]}

EXPATH="${SIMULATIONDIR}/${EXECUTABLE}"

PHI=${parArr[1]}  
cH_2=${parArr[2]}
cP_2=${parArr[3]}
SEED=${parArr[4]}
PREFIX=${parArr[5]}
SUFFIX=${parArr[6]}
INIT=${parArr[7]}



echo "EXECUTABLE:   ${EXECUTABLE} "
echo "EXPATH:   ${EXPATH}"
echo "PHI:  ${PHI}" 
echo "cH_2: ${cH_2}" 
echo "cP_2: ${cP_2}"
echo "SEED: ${SEED}"
echo "PREFIX:   ${PREFIX}"
echo "SUFFIX:   ${SUFFIX}"
echo "INIT: ${INIT}"

# Change to the directory with the compiled source code
cd "${SIMULATIONDIR}"

echo "### Started running simulations."

Rscript --vanilla Run_simulations.R \
    --cH_2 "${cH_2}" \
    --cP_2 "${cP_2}" \
    --executable "${EXPATH}" \
    --phi "${PHI}" \
    --suffix "${SUFFIX}" \
    --parent_outdir "${RESULTSDIR}" \
    --block1 ${BLOCK1} \
    --block2 ${BLOCK2} \
    --init "${INIT}" \
    --sigma ${SIGMA} \
    --betaH "${BETAH}" \
    --betaP "${BETAP}" \
    --prefix "${PREFIX}" \
    --seed "${SEED}" \
    --excheck "${EXCHECK}" \
    --gsize small
```






















