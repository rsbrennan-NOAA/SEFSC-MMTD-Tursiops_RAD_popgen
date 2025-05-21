#!/bin/bash
#SBATCH --job-name=moments_initialRun
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH --mail-type=ARRAY_TASKS
#SBATCH -o %x_%A_%a.out  # %A is the job array ID, %a is the task ID
#SBATCH -e %x_%A_%a.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --ntasks=1 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --partition=standard
#SBATCH --time=5-00
#SBATCH --array=07

source ~/.bashrc
mamba activate moments

# Change to the analysis directory and run the selected model
cd ~/Tursiops-RAD-popgen/analysis/moments/

MODEL_NUMBER=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
REP_NUMBER="assym-2"

python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/moments_model_initialRun.py $MODEL_NUMBER $REP_NUMBER
