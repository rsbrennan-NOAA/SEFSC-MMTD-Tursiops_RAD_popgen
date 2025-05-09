#!/bin/bash
#SBATCH --job-name=m_04-Simple
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --ntasks=1 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --partition=standard
#SBATCH --time=7-00

source ~/.bashrc
mamba activate moments

cd ~/Tursiops-RAD-popgen/analysis/moments/

python -u ~/Tursiops-RAD-popgen/scripts/demography/moments/moments_model_04_Simple.py

