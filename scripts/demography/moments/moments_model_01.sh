#!/bin/bash
#SBATCH --job-name=moments_01
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --partition=standard
#SBATCH --time=72:00:00

source ~/.bashrc
mamba activate moments

cd ~/Tursiops-RAD-popgen/analysis/moments/

python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/moments_model_01.py

