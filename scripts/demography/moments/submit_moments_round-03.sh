#!/bin/bash
#SBATCH --job-name=moments_r3
#SBATCH -o %x_%A_%a.out  # %A is the job array ID, %a is the task ID
#SBATCH -e %x_%A_%a.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --ntasks=1 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --time=3-00
#SBATCH --array=1-192%20  # 8 models × 24 reps = 192 total jobs

source ~/.bashrc
mamba activate moments

# Change to the analysis directory
cd ~/Tursiops-RAD-popgen/analysis/moments/

# Calculate model and rep numbers from array task ID
MODEL_NUMBER=$(( (($SLURM_ARRAY_TASK_ID - 1) % 8) + 1 ))
REP_NUMBER=$(( (($SLURM_ARRAY_TASK_ID - 1) / 8) + 1 ))

MODEL_NUMBER=$(printf "%02d" $MODEL_NUMBER)
REP_NUMBER=$(printf "%02d" $REP_NUMBER)

echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Model number: $MODEL_NUMBER"
echo "Rep number: $REP_NUMBER"

# Run the Round 2 optimization
python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/moments_model_round-03.py $MODEL_NUMBER $REP_NUMBER