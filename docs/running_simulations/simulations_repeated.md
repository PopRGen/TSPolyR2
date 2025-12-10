 ## Pipeline to investigate the effect of initial values

This is the pipeline to investigate the effect of initial genotype frequencies on the dynamics for a fixed combination of parameter values, where only the initial values of genotype frequencies are changed.

The submission script below outlines how $r=50$ repeated simulations for each combination of $\Omega_H=\Omega_M \in \{0.1, 0.2, 0.3\}$, $\Omega_P \in \{0.1, 0.2, 0.3\$, $\Xi_H=\Xi_M=\Xi_P=3$ and $\phi=0.5$ can be run.  A copy of the used simulation schedule can be found [here](../../auxillaries/schedule_ext_reps_CH2_3_CP2_3_phi_0.5_2025-01-07_random.txt).

 Make sure that the submission script below is placed into a directory, which contains a copy of this respository and a `schedules_r50` subdirectory with the simulation schedule. Otherwise the submission script needs to be modified. 
 
If you would like to reproduce the data underlying Figure S6 using the simulation schedule [here](../../auxillaries/schedule_ext_reps_CH2_3_CP2_3_phi_0.5_2025-01-07_random.txt).
 
 
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
 