#!/bin/bash
#SBATCH --job-name=mts_boot
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --partition=standard
#SBATCH --time=7-00

source ~/.bashrc
mamba activate moments

python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/momemts_bootstrap_confInt.py

