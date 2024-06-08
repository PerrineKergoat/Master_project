#!/bin/bash -l
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Isl_rep_analysis
#SBATCH -o out/%x_%A_%a.stdout
#SBATCH -e out/%x_%A_%a.stderr
#SBATCH --mem=25GB
#SBATCH --mail-user perrine.kergoat@unil.ch
#SBATCH --mail-type ALL

module load gcc
module load r

path_data="/work/FAC/FBM/DEE/jgoudet/default/pkergoat/Replicates/Island_model/2pop/data/"

echo "Lancement script 4"
Rscript ./Replicates.R $1 $2 $3 $path_data
echo "Fin script 4"
