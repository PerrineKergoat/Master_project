# Master project: development and application of a novel method to detect selection with genomic data 

This repository contains the scripts and a user guide of the beta method implemented to detect selection from genomic data only. This project was carried out as part of my Master of Science in Behaviour, Evolution and Conservation Computational, Ecology and Evolution project.

The scripts presented are made to simulate a population decomposed in 2, 5 or 10 sub-populations linked by migration and analyse with the beta method the genomes of 1000 individuals among this population. 

## User guide

### Repository architecture

The simulation-beta analysis process in launch using the "Run_script.sh" it requires that all scripts are in a same repository using a SLURM cluster. In addition to the scripts, three directories are required 
- data/
- BED/
- GDS/

A rec_rate_out repository is also needed if the recombination rate analysis is made.

<img width="326" alt="Screenshot 2024-06-08 at 16 53 29" src="https://github.com/PerrineKergoat/Master_project/assets/115696643/a33d30a3-4c33-4347-91e7-a16e85a41197">

### Scripts description

#### Run script

The [Run_script.sh](https://github.com/PerrineKergoat/Master_project/blob/c7b0bb6feeb6b8c0ee1b228030d67950cc566a56/Scripts/Run_script.sh) is the script that coordinates the launch of all the simulation and beta analysis scripts. It must be run with three parameters: 
1. the number of sub-populations simulated
2. the selection coefficient
3. the recombination rate.

Example: 

In the following example we simulate and analyse a population composed of 2 sub-populations with a selection coefficient for mutation under selection of s = 0.030 and a recombination rate of 1.e-8. 

`sbatch Run_script.sh 2 0.030 1e-8`

#### Simulation of populations 

Three SLiM scripts are provided each one corresponding to a different scenario
 
 1- [Sim_2pop_1sel.slim](https://github.com/PerrineKergoat/Master_project/blob/c7b0bb6feeb6b8c0ee1b228030d67950cc566a56/Scripts/Simulations/Sim_2pop_1sel.slim): 

