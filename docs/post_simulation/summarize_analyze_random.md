# summarize_analyze_random.R

## Overview

This script can be used to summarize the maintained genotypes, NLRs, NLR polymorphism for a set of simulations.

## Dependencies

- tidyverse
- argparser
- ggpubr
- grid
- paletteer
- gridExtra

## Usage

```{bash}
Rscript --vanilla summarize_Rgenes.R
```

## Options

* `--prefix`: Constant prefix added to all simulation outputs. Default:  randomcombs. Type: character

* `--suffix`:  suffix for the simulations. Default: 2025-01-27. Type: character
* `--nreps`:  number of simulations. Default: 10000. Type: numeric
* `--parent_dir`:  Directory which has all the results. Default = "/scratch/hm2840/modelling_rerun/results_random". Type: character.
* `--cH_2`:  Value of CH2. Default: 3. Type: numeric
* `--cP_2`:  Value of CP2. Default: 3. Type: numeric
* `--rseed`: Seed. Default: 1700. Type: numeric
* `--phi`: Phi value. Default: 0.5. Type: numeric


## Output

* Note: `<outdir>` is constructed as follows: `<parent_dir>/<prefix>_CH2_<cH_2>_CP2_<cP_2>_random_combis_<suffix>_phi_<phi>`

* Note: `<fileprefix>` is constructed as follows: `CH2_<cH_2>_CP2_<cP_2>_random_r_<nreps>_RSEED_<rseed>`

* `<outdir>/outP_<fileprefix>_phi_<phi>.tsv`: File summarizing the maintained genotypes and alleles for the pathogen for all `<nreps>` simulations. The file has the following columns:

    | Column name | Description |
    | :--- | :--- |
    | phi | Value of phi used for the simulations|
    | CH1 | Value randomly for CH1 for the specific simulation|
    | CP1 |	Randomly drawn value for CP1 for the specific simulation | 
    | seed | seed used for the simulation with the C++ program |
    | Species | Species analyzed |
    | genotypes_remain | summary of genotype which are still present at time 100,000 |
    | range_summary | Numeric value summarizing the maintained ranges. 1: (0), 2: (1), 3: (0+1), 4: (2), 5: (0+2), 6: (1 + 2), 7: (0 + 1 + 2), 8: (3), 9: (0 + 3), 10: (1 + 3), 11: (0 + 1 + 3), 12: (2 + 3), 13: (0 + 2 + 3), 14: (1 + 2 + 3), 15: (0 + 1 + 2 +3) |
    | R_1 | has the virulent allele been maintained at the first locus. TRUE: yes, FALSE: no  |
    | R_2 | has the virulent allele been maintained at the second locus. TRUE: yes, FALSE: no |
    | R_3 |	has the virulent allele been maintained at the third locus. TRUE: yes, FALSE: no|
    | poly_1 | has polymorphism been maintained at the first locus. 0: no, 1: yes |
    | poly_2 | has polymorphism been maintained at the second locus. 0: no, 1: yes  |
    | poly_3 | has polymorphism been maintained at the third locus. 0: no, 1: yes  | 
    | R_alleles | numeric value summarizing which virulence alleles have been maintained. 0: none, 1: only at the first locus, 2: only at the second locus, 3: first and second locus, 4: only third locus, 5: first and third locus, 6: second and third locus, 7: at all loci | 
    | poly | numeric value summarizing at which loci polymorphism has been maintained. 0: none, 1: only at the first locus, 2: only at the second locus, 3: first and second locus, 4: only third locus, 5: first and third locus, 6: second and third locus, 7: at all loci |
    | CH2	| value of $c_H^{(2)}$ used for the simulation |
    | CP2 | value of $c_P^{(2)}$ used for the simulation |

* `<outdir>/outH_<fileprefix>_phi_<phi>.tsv`: File summarizing the maintained genotypes and alleles for both hosts for all `<nreps>` simulations. The file has the following columns:

    | Column name | Description |
    | :--- | :--- |
    | phi | Value of phi used for the simulations|
    | CH1 | Value randomly for CH1 for the specific simulation|
    | CP1 |	Randomly drawn value for CP1 for the specific simulation | 
    | seed | seed used for the simulation with the C++ program |
    | Species | Species analyzed |
    | genotypes_remain | summary of genotype which are still present at time 100,000 |
    | range_summary | Numeric value summarizing the maintained ranges. 1: (0), 2: (1), 3: (0+1), 4: (2), 5: (0+2), 6: (1 + 2), 7: (0 + 1 + 2)|
    | R_1 | has the resistance allele been maintained at the first locus. TRUE: yes, FALSE: no  |
    | R_2 | has the resistance allele been maintained at the second locus. TRUE: yes, FALSE: no |
    | R_3 |	has the resistance allele been maintained at the third locus. TRUE: yes, FALSE: no|
    | poly_1 | has polymorphism been maintained at the first locus. 0: no, 1: yes |
    | poly_2 | has polymorphism been maintained at the second locus. 0: no, 1: yes  |
    | poly_3 | has polymorphism been maintained at the third locus. 0: no, 1: yes  | 
    | R_alleles | numeric value summarizing which resistance alleles have been maintained. 0: none, 1: only private, 2: only shared R allele, 3: private and shared| 
    | poly | numeric value summarizing at which loci polymorphism has been maintained. 0: none, 1: only at the private resistance locus, 2: only shared resistance locus, 3: at the shared and the private resistance locus|
    | CH2	| value of $c_H^{(2)}$ used for the simulation |
    | CP2 | value of $c_P^{(2)}$ used for the simulation |


* `<outdir>/rplt_<fileprefix>_phi_<phi>.pdf`: Scatter plot with two columns (left: host H, right: host M). Each plot depicts the combination of R alleles maintained for each of the `<nreps>` random combinations of $c_H^{(1)}$ (x-axis) and $c_P^{(1)}$ (y-axis). 

* `<outdir>/r_gene_maintained_<fileprefix>_phi_<phi>.pdf`: Plot with two columns (left: host H, right: host M) and two rows (top: maintained R-alleles, bottom: maintained polymorphism). Each plot depicts the results for the `<nreps>` random combinations of $c_H^{(1)}$ (x-axis) and $c_P^{(1)}$ (y-axis).

* `<outdir>/polyplt_<fileprefix>_phi_<phi>.pdf`: Scatter plot with two columns (left: host H, right: host M). Each plot depicts the maintained polymorphism for each of the `<nreps>` random combinations of $c_H^{(1)}$ (x-axis) and $c_P^{(1)}$ (y-axis). 

* `<outdir>/prange_<fileprefix>_phi_<phi>.pdf`: Scatter plot which shows the maintained ranges in the pathogen for each of the the `<nreps>` random combinations of $c_H^{(1)}$ (x-axis) and $c_P^{(1)}$ (y-axis). 

* `<outdir>/ranges_maintained_<fileprefix>_phi_<phi>.pdf`: Composite plot. Top rows shows the maintained resistance allele ranges in host H (left) and host M (right) for all different random combinations of $c_H^{(1)}$ (x-axis) and $c_P^{(1)}$. Bottom plot shows the results of the number of virulence alleles maintained for the pathogen.

* `<outdir>/rbarplt_<fileprefix>_phi_<phi>.pdf`: Bar plot of the proportion of simulations where a given combination of resistance alleles has been maintained in host H (left) and host M (right).

* `<outdir>/polybarplt_<fileprefix>_phi_<phi>.pdf`: Bar plot of the proportion of simulations where a given combination of polymorphism has been maintained in host H (left) and host M (right).

* `<outdir>/bar_r_gene_maintained_<fileprefix>_phi_<phi>.pdf`: Composite plot. Left plot shows the proportion of simulations where a given combination of R-alleles have been maintained for host H (left) and host M (right). The right panel shows the maintained polymorphism for host H (left) and host M (right).

* `<outdir>/rangebarplt_<fileprefix>_phi_<phi>.pdf`: Bar plot of the proportion of simulation where a given combination of ranges (number of R-alleles in the maintained genotypes) has been maintained in host H (left) and host M (right).

* `<outdir>/pbarplt_<fileprefix>_phi_<phi>.pdf`: Bar plot of the proportion of simulation where a given combination of ranges (number of V-alleles in the maintained genotypes) has been maintained in the pathogen. 

* `<outdir>/bar_ranges_maintained_<fileprefix>_phi_<phi>.pdf`: Composite plot. Left plot shows the proportion of simulation where a given combination of ranges (number of resistance alleles) has been maintained in host H (left) and host M (right). The right panel shows the results for the pathogen.

