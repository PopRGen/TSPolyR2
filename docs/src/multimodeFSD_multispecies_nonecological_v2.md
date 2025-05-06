---
title: "README Two host/generalist pathogen model"
author: "Hanna Maerkle"
date: "3/31/2025"
output: html_document
---

<style type="text/css">
  body{
  font-size: 8pt;
}
</style>


# multimodeFSD_multispecies_nonecological_v2.cpp

## Overview

This is the C++-source code for simulating the model dynamics. The source code from Ashby and Boots 2017 for the non-ecological model was used as a starting point. The source code has been extended to accomodate both host species, a command line argument parser has been added and several additional functions have been included to make the code more versatile. Note in the notation of the paper host $H1$ in the code corresponds to host species $H$ and host $H2$ corresponds to host $M$ in the paper

## Dependencies

## Usage

```{bash}
# Compile the source code
make
# Run the program
./Non_ecological_clean --cH_2 3 --cH_1 0.5
```

## Options


| Command line argument | short | Description | Default | Accepted range | Datatype |
| -------| --- | ------------- | --- | -----------| ----- |
| --cP_1 | -p |  Maximum cost of virulence for full range | 0.15 | $0 \leq x \leq 1 $ | double |
| --cP_2 | -q |Parameter defining shape of genotype cost function pathogen ( >0: accelerating, 0: linear, <0: decelerating) | 3 | $-3 \leq x \leq 3 $ | double |
| --cH_1 | -g | Maximum cost of resistance for full range  | 0.15 | $0 \leq x \leq 1 $ | double |
| --cH_2 | -i | Parameter defining the shape of the genotype cost function in both hosts (>0:accelerating, =0:linear, <0:decelerating) | 3 | $-3 \leq x \leq 3 $ | double |
| --sigma | -s | Sigma |  0.85 | $0 \leq x \leq 1 $ | double |
| --betaH | -d | Effect of infection on host  | 1 | $0 \leq x \leq 1 $ | double |
| --betaP | -e | Effect of infection on pathogen | 1 | $0 \leq x \leq 1 $ | double |
| --block1 | -b | null-allele host 1 | 0 | $x \in \{0,1\}$ | int |
| --block2 | -c | null-allele  host 2| 1 | $x \in \{0,1\}$ | int |
| --phi | -a | Relative proportion of host 1| 0.2 | $0 \leq x \leq 1 $ | double | 
| --fprefix | -f | prefix to append to output file | Non_ecological | | string | 
| --seed | -r | seed for simulation. Note choice of seed does not affect simulations for Non_ecological_phi | 1| | int |
| --init | -z | How should genotype frequencies be initialized.| d | 'r', 'f', 'd' |string | 
| --excheck | -y | Should checks be performed if any of the frequencies dropped below 1e-6 and should be reset to 0. 0: do not check, 1: check | 0 | 0 or 1 |  int |
| --help | -h | print the help | | | 

## Output



* The file **"_genotypeH.txt"** file contain the following output:

    | Host | ID | Genotype | P_000 | P_100 | P_010 | P_110 | P_001 | P_101 | P_011 |  P_111 | Basic_fitness | Range | 
    | ---- | ---- | --- | --- | --- | --- | --- | --- | ---  | --- | --- | --- | --- |
    | H1 | 0 | 000 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1  | 1 | 1 | 
    | H1 | 1 | 010 | 0.85 | 0.85 | 1 | 1 | 0.85 | 0.85 | 1 | 1  | 0.883333 | 0.500000 | 
    | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ...  | ... | ... | 
    | H2 | 2 | 001 | 0.85 | 0.85 | 0.85 | 0.85 | 1 | 1 | 1 | 1  | 0.883333 | 0.500000 |
    | H2 | 3 | 101 | 0.7225 | 0.85 | 0.7225 | 0.85 | 0.85 | 1 | 0.85 | 1  | 0.76667 | 0.500000 |

    The output should be read as follows:
    - Host: Host species 
    - ID: ID of the genotype in the host species given in column Host
    - Genotype: Genotype of the host
    - P_000: Chance of infection when interacting with pathogen genotype P_000
    - P_100: Probability of infection when interacting with pathogen genotype P_100
    - ...
    - Basic_fitness: Fitness in the absence of the pathogen
    - Range: Resistance range of the host genotype

* The file **"_genotypeP.txt"** contains the following output:

    | Para | ID | Genotype | H1_000 | H1_010 | H1_001 | H1_011 | H2_000 | H1_100 | H2_001 |  H2_101 | Basic_fitness | Range | 
    | ---- | ---- | --- | --- | --- | --- | --- | --- | ---  | --- | --- | --- | --- |
    | P | 0 | 000 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1  | 1 | 1 | 
    | P | 1 | 100 | 0.85 | 0.85 | 1 | 1 | 0.85 | 0.85 | 1 | 1  | 0.883333 | 0.500000 | 
    | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ...  | ... | ... | 
    | P | 6 | 011 | 1 | 1 | 1 | 1 | 1 | 0.85 | 1 | 0.85 | 0.966524 | 0.666667 |
    | P |7  | 111 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0.9 | 1.000000 |

    The output should be read as follows:
        
    - Host: Host species 
    - ID: ID of the genotype in the host species given in column Host
     - Genotype: Genotype of the host
    - P_000: Chance of infection when interacting with pathogen genotype P_000
    - P_100: Probability of infection when interacting with pathogen genotype P_100
    - ...
    - Basic_fitness: Fitness in the absence of the pathogen
    - Range: Resistance range of the host genotype

* The file **"_dynamics.txt"** contains all frequency dynamics in the simulation.



## Additional information 

### Solving of the differential equations

The differential equations are solved using the Cash-Karp method (Runge-Kutta method with adaptive step size) which estimates the forth-order and five-order accurate solution. It defines a new step size based on estimating the error of the forth-order accurate solution by comparing it to the fifth-order accurate solution. The original code in Ashby and Boots 2017 followed the numerical recipes on implementing the Cash-Karp method with adaptive step size.

The Butcher-tableau used looks as follows:

$$
\begin{matrix} 
    & a_i & & & b_{ij} & & & c_i & c^{\star}_i \\
i=1 & 0 & & & & & & \frac{37}{378} & \frac{2825}{27684}\\
i=2 & 1/5 & \frac{1}{5} & & & &  & 0 & 0\\
i=3 & 3/10 & \frac{3}{40} & \frac{9}{40} & & &  & \frac{250}{621} & \frac{18575}{48384} \\
i=4 & 3/5 & \frac{3}{10} & \frac{-9}{10} & \frac{6}{5} & & & \frac{125}{594} & \frac{13525}{55296} \\
i=5 & 1 & -\frac{11}{54} & \frac{5}{2} & -\frac{70}{27} &  \frac{35}{27} &  & 0 & \frac{277}{14336} \\
i=6 & 7/8 & \frac{1631}{55296} & \frac{175}{512} & \frac{575}{13824} & \frac{44275}{110592} &  \frac{253}{4096} & \frac{512}{1771} & \frac{1}{4} \\
    &     & j=1 & j=2 & j=3 & j=4 & j=5 & &
\end{matrix}
$$

This gives the following coefficients for calculating the errors:

$dc_1 = c_1 - c^{\star}_1 = c_1 - \frac{2825.0}{27648.0}$  
$dc_2 = c2 - c^{\star}_2= 0$  
$dc_3 = c_3 - c^{\star}_3 = c_3 - \frac{18575}{48384.0}$  
$dc_4 = c_4 - c^{\star}_4 = c_4 - \frac{13525}{55296.0}$  
$dc_5 = c_5 - c^{\star}_5 = -\frac{277.00}{14336}$  
$dc_6 = c_6 - c^{\star}_6 = c_6 - 0.25$

The error for the j-th ODE is calcuated as:

$\Delta_{1,j} = h * \sum_{i=1}^6 dc_i k_{i,j}$

The desired accuracy for the j-th ODE is calculated as:

$\Delta_{0,j} = \epsilon \cdot (|y_j |  + | h * k_{1,j} | + \text{TINY})$ where $y_j$ is the current value of the j-th dependent variable and $k_{1,j}$ is derivative of the k-th variable with respect to time.

For all simulations $\epsilon = 1e^{-6}$ and $\text{TINY} = 1e^-6$

The stepsize is adjusted as follows:

$$
h_\text{next} = 
\begin{cases}
S * h * |\frac{\Delta_1}{\Delta_0}|^{-0.2} & \quad \text{when $\frac{\Delta_1}{\Delta_0} \leq 1 $}  \\ 
S * h * |\frac{\Delta_1}{\Delta_0}|^{-0.25} & \quad \text{when $\frac{\Delta_1}{\Delta_0} > 1 $}
\end{cases}
$$

Note if the stepsize is decreased, the stepsize will not be decreased by more than a factor of 10. Note is $\frac{\Delta_1}{\Delta_0} < {\frac{5}{S}}^{\frac{1}{-0.2}}$ then the new step size is calculated as: $h_\text{next} = 5 * h$
