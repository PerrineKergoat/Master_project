#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=W500snps
#SBATCH -o out/%x_%j.stdout
#SBATCH -e out/%x_%j.stderr
#SBATCH --mem=6GB
#SBATCH --mail-user perrine.kergoat@unil.ch
#SBATCH --mail-type ALL

module load gcc
module load slim/4.0.1
module load python
module load r
source /work/FAC/FBM/DEE/jgoudet/default/pkergoat/pyslim_venv/bin/activate

echo "Lancement script 1"
slim -d replicate_number=$SLURM_JOB_ID Sim_2pop_1sel.slim
echo "Fin script 1"

echo "Lancement script 2"
python ./Recapitation.py $SLURM_JOB_ID
echo "Fin script 2"

echo "Lancement script 3"
Rscript ./Comparison.R $SLURM_JOB_ID
echo "Fin script 3"

echo "Lancement script 4"
Rscript ./FST_analysis.R $SLURM_JOB_ID
echo "Fin script 4"
