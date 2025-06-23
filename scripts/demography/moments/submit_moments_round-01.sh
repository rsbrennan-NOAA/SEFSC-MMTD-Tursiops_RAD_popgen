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
#SBATCH --time=4-00
#SBATCH --array=1-192%18  # 8 models × 12 reps = 96 total jobs

source ~/.bashrc
mamba activate moments

# Change to the analysis directory and run the selected model
cd ~/Tursiops-RAD-popgen/analysis/moments/

# Calculate model and rep numbers from array task ID
MODEL_NUMBER=$(( (($SLURM_ARRAY_TASK_ID - 1) % 8) + 1 ))
REP_NUMBER=$(( (($SLURM_ARRAY_TASK_ID - 1) / 8) + 1 ))

MODEL_NUMBER=$(printf "%02d" $MODEL_NUMBER)
REP_NUMBER=$(printf "%02d" $REP_NUMBER)

echo "model number: $MODEL_NUMBER"
echo "rep number: $REP_NUMBER"

#python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/moments_model_initialRun.py $MODEL_NUMBER $REP_NUMBER
python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/moments_model_round-01_m${MODEL_NUMBER}.py $MODEL_NUMBER $REP_NUMBER
