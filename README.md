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
