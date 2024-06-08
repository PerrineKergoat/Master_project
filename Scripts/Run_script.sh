#!/bin/bash -l
#SBATCH --time=15:00:00
#SBATCH --array=1-120
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Isl_mod
#SBATCH -o out/%x_%A_%a.stdout
#SBATCH -e out/%x_%A_%a.stderr
#SBATCH --mem=3GB
#SBATCH --mail-user perrine.kergoat@unil.ch
#SBATCH --mail-type ALL

module purge
dcsrsoft use arolle

module load gcc
module load slim/4.0.1
module load python
module load r
source /work/FAC/FBM/DEE/jgoudet/default/pkergoat/pyslim_venv/bin/activate

gen_600="600"
gen_800="800"
gen_1000="1000"
path_data="/work/FAC/FBM/DEE/jgoudet/default/pkergoat/Replicates/Island_model/2pop/data/"

echo "Lancement script 1: SLiM"
slim -d job_id=$SLURM_ARRAY_JOB_ID -d rep_nb=$SLURM_ARRAY_TASK_ID -d nb_pop=$1 -d sel_coeff=$2 -d recomb_rate=$3 -d "path='$path_data'" Sim_2pop_1sel.slim
echo "Fin script 1"

echo -e "\nLancement script 2: Recapitation.py for 600 generations"
python ./Recapitation.py $SLURM_ARRAY_JOB_ID  $SLURM_ARRAY_TASK_ID $3 $gen_600 $1 $path_data
echo "Fin script 2 600 generations"

echo -e "\nLancement script 2: Recapitation.py for 800 generations"
python ./Recapitation.py $SLURM_ARRAY_JOB_ID  $SLURM_ARRAY_TASK_ID $3 $gen_800 $1 $path_data
echo "Fin script 2_800"

echo -e "\nLancement script 2: Recapitation.py for 1000 generations"
python ./Recapitation.py $SLURM_ARRAY_JOB_ID  $SLURM_ARRAY_TASK_ID $3 $gen_1000 $1 $path_data
echo "Fin script 2_1000"

echo -e "\nLancement script 3: Comparison.R for 600 generations"
Rscript ./Comparison.R $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID $1 $gen_600 $path_data
echo "Fin script 3_600"

echo -e "\nLancement script 3_800"
Rscript ./Comparison.R $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID $1 $gen_800 $path_data
echo "Fin script 3_800"

echo -e "\nLancement script 3_1000"
Rscript ./Comparison.R $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID $1 $gen_1000 $path_data
echo "Fin script 3_1000"
