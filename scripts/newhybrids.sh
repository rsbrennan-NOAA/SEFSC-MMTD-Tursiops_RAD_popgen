#!/bin/bash
#SBATCH --job-name=newhybrids
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -n 1
#SBATCH --mem=40G
#SBATCH --partition=standard
#SBATCH --time=100:00:00
#SBATCH --array=1-10


cd ~/Tursiops-RAD-popgen/analysis/pop_structure/newhybrids

mkdir -p ~/Tursiops-RAD-popgen/analysis/pop_structure/newhybrids/run_${SLURM_ARRAY_TASK_ID}

cd ~/Tursiops-RAD-popgen/analysis/pop_structure/newhybrids/run_${SLURM_ARRAY_TASK_ID}

echo "Job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} started at: $(date)"

# Run newhybrids
~/bin/newhybrids/newhybrids-no-gui-linux.exe --no-gui \
  --data-file ~/Tursiops-RAD-popgen/analysis/variants/all_variants.3.ped \
  --burn-in 100000 \
  --num-sweeps 1000000 

cp aa-PofZ.txt ~/Tursiops-RAD-popgen/analysis/pop_structure/newhybrids/aa-PofZ_${SLURM_ARRAY_TASK_ID}.txt

echo "Job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} completed at: $(date)"
