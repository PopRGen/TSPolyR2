## Recipe 2: Running simulations for random combinations of $\Omega_P$, $\Omega_H$ for fixed values of $\Xi_H$, $\Xi_P$ and $\phi$.


Use the script [src/Run_simulation_randomdraw.R](./docs/src/Run_simulation_randomdraw.md) to run a fixed number of simulations drawing random values for $\Omega_H$, $\Omega_P$ and random initial genotype frequencies. Note previously the paramters were named as follows:

| New name | Old name | Description |
| :---- | :----| :--- |
| $\Omega_H$  | $c_H^{(1)}$ | maximum cost of resistance |
| $\Omega_P$  | $c_P^{(1)}$ | maximum cost of virulence |
| $\Xi_H$  | $c_H^{(2)}$ | parameter for shape of resistance cost function |
| $\Xi_P$  | $c_P^{(2)}$ | parameter for shape of virulence cost function |

### Example 1 running 100 simulations with random initial genotype frequencies and randomly drawing values for $\Omega_H=\Omega_M$ and $\Omega_P$.

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
This runs 100 simulations. For each simulation $\Xi_H = \Xi_P = 3$ and $\phi=0.5$. The values of $\Omega_H$ for each simulation is drawn from $\mathcal{U}(0.01,0.3)$ and the value of $\Omega_P$ for each simulation is drawn from $\mathcal{U}(0.01,0.3)$.

## Running a subset of simulations underlying Figures 2+3+4, S4, S5 and S7

The job submission script below can be used to generate 10,000 random simulations for a fixed combination of shape parameters $\Xi_H = \Xi_M = \Xi_P = 3$, using a uniform prior $\mathcal{U}(0.01,0.3)$ for $\Omega_H=\Omega_M$ and $\mathcal{U}(0.01,0.3)$ for $\Omega_P$ with random initial genotype frequencies. The script will generate 10,000 simulations for different values of the relative proportion of host H ($\phi \in \{0.5,0.7,0.9,1.0\}$).

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
