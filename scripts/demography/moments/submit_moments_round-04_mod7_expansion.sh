#!/bin/bash
#SBATCH --job-name=moments_inter
#SBATCH -o %x_%A_%a.out  # %A is the job array ID, %a is the task ID
#SBATCH -e %x_%A_%a.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --ntasks=1 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --time=20:00:00
#SBATCH --array=1-20%20  # 8 models × 10 reps = 80 total jobs

source ~/.bashrc
mamba activate moments

# set parameter here. 
ROUND_NUMBER=4
PERTURB=1
MAXITER=5
METHOD="fmin"
SFS_FILE="/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_equalSS/dadi/Coastal_Atlantic-Coastal_Gulf-Intermediate-Offshore.sfs"

# Change to the analysis directory
cd ~/Tursiops-RAD-popgen/analysis/moments/

# Calculate model and rep numbers from array task ID
MODEL_NUMBER="07"
REP_NUMBER=$SLURM_ARRAY_TASK_ID

MODEL_NUMBER=$(printf "%02d" $MODEL_NUMBER)
REP_NUMBER=$(printf "%02d" $REP_NUMBER)

echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Model number: $MODEL_NUMBER"
echo "Rep number: $REP_NUMBER"
echo "Round: $ROUND_NUMBER"
echo "SFS file: $SFS_FILE"
echo "Parameters: PERTURB=$PERTURB, MAXITER=$MAXITER, METHOD=$METHOD"

# Run the optimization
python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/moments_mod7_expansion_round-03_plus.py $MODEL_NUMBER $REP_NUMBER $ROUND_NUMBER $PERTURB $MAXITER $METHOD $SFS_FILE