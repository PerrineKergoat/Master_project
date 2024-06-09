# Master project: development and application of a novel method to detect selection with genomic data 

This repository contains the scripts and a user guide of the beta method implemented to detect selection from genomic data only. This project was carried out as part of my Master of Science in Behaviour, Evolution and Conservation Computational, Ecology and Evolution project.

The scripts presented are made to simulate a population decomposed in 2, 5 or 10 sub-populations linked by migration and analyse with the beta method the genomes of 1000 individuals among this population. 

# User guide

## Repository architecture

The simulation-beta analysis process in launch using the "Run_script.sh" it requires that all scripts are in a same repository using a SLURM cluster. In addition to the scripts, three directories are required 
- data/
- BED/
- GDS/

A rec_rate_out repository is also needed if the recombination rate analysis is made.

<img width="326" alt="Screenshot 2024-06-08 at 16 53 29" src="https://github.com/PerrineKergoat/Master_project/assets/115696643/a33d30a3-4c33-4347-91e7-a16e85a41197">

## Scripts description

### Run script

The [Run_script.sh](https://github.com/PerrineKergoat/Master_project/blob/c7b0bb6feeb6b8c0ee1b228030d67950cc566a56/Scripts/Run_script.sh) is the script that coordinates the launch of all the simulation and beta analysis scripts. It must be run with three parameters: 
1. the number of sub-populations simulated
2. the selection coefficient
3. the recombination rate.

Example: 

In the following example we simulate and analyse a population composed of 2 sub-populations with a selection coefficient for mutation under selection of s = 0.030 and a recombination rate of 1.e-8. 

`$ sbatch Run_script.sh 2 0.030 1e-8`

### Simulation of populations with SLiM

Three SLiM scripts are provided each one corresponding to different scenario of simulations.
 
 1. [Sim_2pop_1sel.slim](https://github.com/PerrineKergoat/Master_project/blob/c7b0bb6feeb6b8c0ee1b228030d67950cc566a56/Scripts/Simulations/Sim_2pop_1sel.slim): this script simulates as many sub-populations as required by the user connecting with migration accroding an island model. In this script, the recombination rate is constant along the genome using the rate givven by the user. 
 2. [Sim_2pop_1sel_recomb.slim](https://github.com/PerrineKergoat/Master_project/blob/c7b0bb6feeb6b8c0ee1b228030d67950cc566a56/Scripts/Simulations/Sim_2pop_1sel_recomb.slim): this script simulates as many sub-populations as required by the user connecting with migration accroding an island model. In this script, the recombination rate is varying, the genome is cut in 100 chunks with recombination rates determined by a gamma distribution.  
 3. [Sim_stepstone.slim](https://github.com/PerrineKergoat/Master_project/blob/c7b0bb6feeb6b8c0ee1b228030d67950cc566a56/Scripts/Simulations/Sim_stepstone.slim): this script simulates as many sub-populations as required by the user connecting with migration accroding an steppingstone model. The population under selection must be manually changed in the script and is not used as a parameter to launch the Run_script.sh. 

The output of these scipt is a `.trees` file stored in the data repository.

### Recapitation of SLiM simulations using Python

The [Recapitation.py](https://github.com/PerrineKergoat/Master_project/blob/ca293c6686c73951234345f36ced3a1fce81b6ea/Scripts/Simulations/Recapitation.py) script creates the evolutionary of the SLiM simulations, it adds the neutral mutation background to the genomes.

The script use the `.trees` files previoously produced and creates a `.vcf` file as an output.

### Beta analysis of a population

The beta analysis of the genomes to detect selection is implemented in the [Comparison.R](https://github.com/PerrineKergoat/Master_project/blob/ca293c6686c73951234345f36ced3a1fce81b6ea/Scripts/beta_analysis/Comparison.R) script. 

This analysis is made based on a `.vcf` file containing all the individuals studied. It produces an `.RData` environment containing all data for the graphical analysis. 

### Analysis of replicates

The different simulations are replicated to ensure reliable results for the Master Project. The analysis of those replicates is launched with the [Summary_replicates_script.sh](https://github.com/PerrineKergoat/Master_project/blob/9e54802e65dfe96ce27c7a8ae53740ac163bad77/Scripts/Replicates_analysis/Summary_replicates_script.sh) script.

Depending if the an analysis of the recombination rate around the beneficial mutation is expected, two scripts can be run with the bash script: 
1. [Replicates.R](https://github.com/PerrineKergoat/Master_project/blob/9e54802e65dfe96ce27c7a8ae53740ac163bad77/Scripts/Replicates_analysis/Replicates.R) that is usefull when no analysis of the recombination rate is desired
2. [Replicates_rec_rate.R](https://github.com/PerrineKergoat/Master_project/blob/9e54802e65dfe96ce27c7a8ae53740ac163bad77/Scripts/Replicates_analysis/Replicates_rec_rate.R) if the supplementary analysis is wanted

This analysis of all the `.RData` files created for the replicated simulation and group them into one `.RData` file.
