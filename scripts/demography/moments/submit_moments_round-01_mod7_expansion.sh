#!/bin/bash
#SBATCH --job-name=moments_initialRun
#SBATCH -o %x_%A_%a.out  # %A is the job array ID, %a is the task ID
#SBATCH -e %x_%A_%a.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --ntasks=1 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --partition=standard
#SBATCH --time=1-00
#SBATCH --array=41-80%20  # 8 models × 12 reps = 96 total jobs

source ~/.bashrc
mamba activate moments

# Change to the analysis directory and run the selected model
cd ~/Tursiops-RAD-popgen/analysis/moments/

# Calculate model and rep numbers from array task ID
REP_NUMBER=$SLURM_ARRAY_TASK_ID

REP_NUMBER=$(printf "%02d" $REP_NUMBER)
MOD_NUMBER="07"

echo "rep number: $REP_NUMBER"

#python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/moments_model_initialRun.py $MODEL_NUMBER $REP_NUMBER
python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/moments_mod7_expansion_round-01.py $MOD_NUMBER $REP_NUMBER