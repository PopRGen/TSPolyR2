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

Find a job submission script below which can be used to generate 10,000 random simulations for $c_H^{(2)}=c_H^{(2)}=3$ and $\phi=0.1$, using a uniform prior $\mathcal{U}(0.01,0.3)$ for $c_P^{(1)}$ and $\mathcal{U}(0.01,0.3)$ for $c_H^{(1)}$ with random initial frequencies.

```{bash,eval=F}
#!/bin/bash

#SBATCH --mem=4GB
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --job-name=run_randomcombs
#SBATCH --output=logs/run_randomcombs_CH2_3_CP2_3_new_%A_%a.out
#SBATCH --error=logs/run_randomcombs_CH2_3_CP2_3_new_%A_%a.err
#SBATCH --tasks=1
#SBATCH --array=1

module load gcc/10.2.0
module load r/gcc/4.4.0

PHIVALUES=(0.5 0.7 0.9 1)
SEEDVALUES=(1700 2700 3700 4700)
ID=${SLURM_ARRAY_TASK_ID}
PHI=${PHIVALUES[ID]}
PHIMIN=${PHI}
PHIMAX=${PHI}
RSEED=${SEEDVALUES[ID]}

SUFFIX="2025-01-27_phi_${PHIMIN}"
PREFIX="randomcombs"
CH2=3
CP2=3
NREPS=10000

# Set the desired parameters for the simulation
BETAP=1
BETAH=1
SIGMA=0.85
BLOCK1=0
BLOCK2=1
EXCHECK=1
CH1MIN=0.01
CH1MAX=0.3
CP1MIN=0.01
CP1MAX=0.3
INIT="r"
NLOCI=3
EXCHECK=1

# These two paths are relative to the working_dir
WORKINGDIR=`pwd`
OUTDIR="${WORKINGDIR}/results_random"
SIMULATIONDIR="${WORKINGDIR}/TSPolyR2/src"
EXECUTABLE="Non_ecological_clean"
EXPATH="${SIMULATIONDIR}/${EXECUTABLE}"
TAR="${OUTDIR}/${PREFIX}_CH2_${CH2}_CP2_${CP2}_random_combis_${SUFFIX}/simulations/CH2_${CH2}_CP2_${CP2}_random_r_${NREPS}_RSEED_${RSEED}.tar"



echo "WORKINGDIR:   ${WORKINGDIR}"
echo "PREFIX:   ${PREFIX}"
echo "SUFFIX:   ${SUFFIX}"
echo "BETAP:    ${BETAP}" 
echo "BETAH:    ${BETAH}" 
echo "SIGMA:    ${SIGMA}"
echo "BLOCK1:   ${BLOCK1}"
echo "BLOCK2:   ${BLOCK2}"
echo "OUTDIR:   ${OUTDIR}"
echo "SIMULATIONDIR: ${SIMULATIONDIR}"
echo "OUTDIR: ${OUTDIR}"
echo "EXECUTABLE: ${EXECUTABLE}"
echo "EXPATH:   ${EXPATH}"
echo "CH1MIN:  ${CH1MIN}" 
echo "CH1MAX:  ${CH1MAX}" 
echo "CP1MIN:  ${CP1MIN}" 
echo "CP1MAX:  ${CP1MAX}" 
echo "PHIMIN:  ${PHIMIN}" 
echo "PHIMAX:  ${PHIMAX}" 
echo "cH_2: ${cH_2}" 
echo "cP_2: ${cP_2}"
echo "RSEED: ${RSEED}"
echo "EXCHECK: ${EXCHECK}"
echo "INIT: ${INIT}"
echo "NLOCI: ${NLOCI}"
echo "NREPS: ${NREPS}"
echo "TAR: ${TAR}"

# Change to the directory with the compiled source code
cd "${SIMULATIONDIR}"

echo "### Started running simulations."

Rscript --vanilla Run_simulation_randomdraw.R \
    --cH_2 "${CH2}" \
    --cP_2 "${CP2}" \
    --executable "${EXPATH}" \
    --cH1_min "${CH1MIN}" \
    --cH1_max "${CH1MAX}" \
    --cP1_min "${CP1MIN}" \
    --cP1_max "${CP1MAX}" \
    --phi_min "${PHIMIN}" \
    --phi_max "${PHIMAX}" \
    --suffix "${SUFFIX}" \
    --prefix "${PREFIX}" \
    --outdir "${OUTDIR}" \
    --block1 ${BLOCK1} \
    --block2 ${BLOCK2} \
    --init "${INIT}" \
    --sigma ${SIGMA} \
    --betaH "${BETAH}" \
    --betaP "${BETAP}" \
    --rseed "${RSEED}" \
    --excheck "${EXCHECK}" \
    --nreps "${NREPS}" \
    --nloci "${NLOCI}" 

# Check if a compressed archive already exists
if [ -f "${TAR}.gz" ]
then
	mv "${TAR}.gz" "${TAR%tar.gz}_old.tar.gz"
fi

# Compress the uncompressed tar archive
pigz -k -p ${SLURM_CPUS_PER_TASK} ${TAR}
# Remove the tar archive
rm ${TAR}
```

### Summarize the results

The results can be summarize with 

This submission script below can be used to summarize the simulations for $c_H^{(2)}=c_P^{(2)}=3$, $\phi=0.5$ and RSEED=1700 (the seed making sure that the random draws in R are reproducible).

```{bash,eval=F}
#!/bin/bash

#SBATCH --time=40:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --error=logs/plot_results_random_CH2_CP2_%A_%a.err
#SBATCH --output=logs/plot_results_random_CH2_CP2_%A_%a.out
#SBATCH --array=1

module load r/gcc/4.4.0

PARENT_DIR="results_random"
PREFIX="randomcombs"
SUFFIX="2025-01-27"
NREPS=10000

CH2VALUES=(3)
CP2VALUES=(3)
PHIVALUES=(0.5)
RSEEDVALUES=(1700)

ID="${SLURM_ARRAY_TASK_ID}"
PHI=${PHIVALUES[ID]}
RSEED=${RSEEDVALUES[ID]}
CH2=${CH2VALUES[ID]}
CP2=${CP2VALUES[ID]}

echo "PARENTDIR:    ${PARENT_DIR}"
echo "PREFIX:       ${PREFIX}"
echo "SUFFIX:       ${SUFFIX}"
echo "NREPS:        ${NREPS}"
echo "ID:           ${ID}"
echo "CH2:          ${CH2}"
echo "CP2:          ${CP2}"
echo "PHI:          ${PHI}"
echo "RSEED:        ${RSEED}"

Rscript --vanilla TSPolyR2/post_simulation/summarize_analyze_random.R \
    --parent_dir ${PARENT_DIR} \
    --prefix ${PREFIX} \
    --suffix ${SUFFIX} \
    --nreps ${NREPS} \
    --cH_2 ${CH2} \
    --cP_2 ${CP2} \
    --phi ${PHI} \
    --rseed ${RSEED}
```

### Plot Figures 2-6

Figures 2-6 can be plotted using the script [post_simulation/plot_figures_2-6.R](./docs/post_simulation/plot_figures_2-6.R). It requires that all pairwise combinations of $c_H^{(2)} \in \{-3,0,3\}$, $c_P^{(2)} \in \{-3,0,3\}$ and $\phi \in \{0.5, 0.7, 0.9, 1\}$ prior to runnint the script. Find a job submission script below.

```{bash, eval=F}
#!/bin/bash

#SBATCH --time=40:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --output=logs/generate_plts_2-6_%A.out
#SBATCH --error=logs/generate_plts_2-6_%A.err

module load r/gcc/4.4.0

Rscript --vanilla TSPolyR2/post_simulation/plot_figures_2-6.R \
	--parentdir `pwd`/results_random \
	--suffix 2025-01-25 \
	--outdir Figures  
```


## Recipe 3: Running simulations for fixed combinations of $c_P^{(1)}$, $c_H^{(1)}$, $c_H^{(2)}$, $c_P^{(2)}$ and $\phi$ with random intial frequencies.

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
The simulation schedules can be generated with the script [auxillaries/write_schedules_r50.R](./docs/auxillaries/write_schedules_r50.R).

### Collect the results across seeds

The script `TSPolyR/post_simulation/investigate_initial_condition.R` can be used for comparing and generating overview of the effect of initial genotype frequency for a fixed parameter combination. Find a job submission script below how to generate these overviews for all pairwise combinations of $c_H^{(1)} \in \{0.1, 0.2, 0.3\}$, $c_P^{(1)} \in \{0.1, 0.2, 0.3\}$, for $c_H^{(2)}=c_P^{(2)}=0$ and $\phi=0.5$ and for all pairwise combinations of $c_H^{(1)} \in \{0.1, 0.2, 0.3\}$, $c_P^{(1)} \in \{0.1, 0.2, 0.3\}$, for $c_H^{(2)}=c_P^{(2)}=3$ and $\phi=0.5$.

```{bash,eval=F}
#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --mem=16GB
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --error=logs/run_seedcomp_phi_0.5_r50_%A_%a.err
#SBATCH --output=logs/run_seedcomp_phi_0.5_r50_%A_%a.out
#SBATCH --array=0-1

module load r/gcc/4.4.0


PHI=0.5
PREFIX="ext_reps"
SUFFIX="2025-01-07"
ID=${SLURM_ARRAY_TASK_ID}

CH2VALUES=(0 3)
CP2VALUES=(0 3)
CH2=${CH2VALUES[ID]}
CP2=${CP2VALUES[ID]}

INDIR="results_r50/${PREFIX}_CH2_${CH2}_CP2_${CP2}_${SUFFIX}/all_dynamics_long"
PLOTDIRTOP="Figures_seed_comparison/${PREFIX}_CH2_${CH2}_CP2_${CP2}_PHI_${PHI}_${SUFFIX}"
NREP=50
SEED_START=400
SEED_INTERVAL=400


CH1VALUES=(0.1 0.2 0.3)
CP1VALUES=(0.1 0.2 0.3)


echo "WDIR:  ${WDIR}"

for j in "${!CH1VALUES[@]}"
do
        CH1=${CH1VALUES[j]}
        for k in "${!CP1VALUES[@]}"
        do
                CP1=${CP1VALUES[k]}
		        PLOTDIR="${PLOTDIRTOP}/CH1_${CH1}_CP1_${CP1}"
                echo "### Generating summaries for:"
                echo "PREFIX:   ${PREFIX}"
                echo "SUFFIX:   ${SUFFIX}"
                echo "INDIR:    ${INDIR}"
                echo "PLOTDIR:  ${PLOTDIR}"
                echo "CH2:      ${CH2}"
                echo "CP2:      ${CP2}"
                echo "CH1:      ${CH1}"
                echo "CP1:      ${CP1}"
		        echo "PHI:      ${PHI}"
                echo "NREP:	${NREP}"
		        echo "SEED_START ${SEED_START}"
		        echo "SEED_INTERVAL ${SEED_INTERVAL}"		

                Rscript --vanilla TSPolyR/post_simulation/investigate_initial_condition.R \
                        --prefix ${PREFIX} \
                        --suffix ${SUFFIX} \
                        --indir ${INDIR} \
                        --plotdir ${PLOTDIR} \
                        --tselect 40000 \
                        --cH_2 ${CH2} \
                        --cP_2 ${CP2} \
                        --cH_1 ${CH1} \
                        --cP_1 ${CP1} \
                        --phi ${PHI} \
			            --nrep ${NREP} \
			            --seed_start ${SEED_START} \
			            --seed_interval ${SEED_INTERVAL}
        done
done 

echo "Compressing the "
tar --use-compress-program="pigz -k -p $SLURM_CPUS_PER_TASK" -cvf ${INDIR}_phi_${PHI}.tar.gz ${INDIR}/*_phi_${PHI}*

```

### Plot dynamics for all simulations with seed=1600

The script `TSPolyR/post_simulation/Plot_fluctuations.R` can be used to plot the dynamics for all simulations where the seed was 1600 for the C++ program. This script was used to generate Figure 1. Find a job submission script below on how to run the script for all pairwise combinations of $c_H^{(1)} \in \{0.1, 0.2, 0.3\}$, $c_P^{(1)} \in \{0.1, 0.2, 0.3\}$, for $c_H^{(2)}=c_P^{(2)}=0$ and $\phi=0.5$ and for all pairwise combinations of $c_H^{(1)} \in \{0.1, 0.2, 0.3\}$, $c_P^{(1)} \in \{0.1, 0.2, 0.3\}$, for $c_H^{(2)}=c_P^{(2)}=3$ and $\phi=0.5$.

```{bash,eval=F}
#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --error=logs/dynamicsplt_%A_%a.err
#SBATCH --output=logs/dynamicsplt_%A_%a.out
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --array=0-1

module load r/gcc/4.4.0

# Script to plot dynamics for seed 1600


ID=${SLURM_ARRAY_TASK_ID}
CH2VALUES=(0 3)
CP2VALUES=(0 3)

CH2=${CH2VALUES[ID]}
CP2=${CP2VALUES[ID]}

SUFFIX="2025-01-07"
PARENTOUTDIR="results_r50"
SEED=1600
PLOTDIR="Figures_seed_comparison/Dynamics_seed_1600"
PREFIX="ext_reps"

echo  "CH2:		${CH2}"
echo  "CP2:		${CP2}"
echo  "PHI:		${PHI}"
echo  "SUFFIX:		${SUFFIX}"
echo  "SEED:		${SEED}"
echo  "PLOTDIR: 	${PLOTDIR}"
echo  "PREFIX:  	${PREFIX}"
echo  "PARENTOUTDIR:	${PARENTOUTDIR}"

PHIVALUES=(0.5)
CH1VALUES=(0.1 0.2 0.3)
CP1VALUES=(0.1 0.2 0.3)

for i in "${!PHIVALUES[@]}"
do
	PHI=${PHIVALUES[i]}
    for j in "${!CH1VALUES[@]}"
    do
		CH1=${CH1VALUES[j]}
		for k in "${!CP1VALUES[@]}"
		do
			CP1=${CP1VALUES[k]}
     			echo "PHI:	${PHI}"
 			echo "CH1:	${CH1}"
			echo "CP1:	${CP1}"
			Rscript --vanilla TSPolyR/post_simulation/Plot_fluctuations.R \
				--phi ${PHI} \
				--suffix ${SUFFIX} \
				--prefix ${PREFIX} \
				--cH_2 ${CH2} \
				--cP_2 ${CP2} \
				--cH_1 ${CH1} \
				--cP_1 ${CP1} \
				--plotdir ${PLOTDIR} \
				--parent_outdir ${PARENTOUTDIR} \
				--seed ${SEED}

		done
    done
done
```

### Plot supplementary Figure 1

Supplementary figure 1 can be generated by running script [post_simulation/plot_supplementary_1.R])(./docs/post_simulation/plot_supplementary_1.md) with default settings. Make sure that you R working directory is set to the path parent directory containing the results of the standard workflow. 




















