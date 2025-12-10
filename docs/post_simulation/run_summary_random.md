## Summarizing simulations for fixed shaped parameters

The following script can be used to summarize a set of x simulations for fixed values of the shape parameters $\Xi_H = \Xi_M$, $\Xi_P$ and a fixed value of $\phi$ with the values of $\Omega_H=\Omega_M$, $\Omega_P$ drawn at random. For details see the manuscript. Note previously the paramters were named as follows:

| New name | Old name | Description |
| :---- | :----| :--- |
| $\Omega_H$  | $c_H^{(1)}$ | maximum cost of resistance |
| $\Omega_P$  | $c_P^{(1)}$ | maximum cost of virulence |
| $\Xi_H$  | $c_H^{(2)}$ | parameter for shape of resistance cost function |
| $\Xi_P$  | $c_P^{(2)}$ | parameter for shape of virulence cost function |

The script below can be extended by adjusting the arrays correspondingly.

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