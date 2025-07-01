#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --partition standard

source ~/.bashrc
conda activate snakemake-9.6.2

cd /home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments
snakemake --profile slurm -s moments.smk --unlock

snakemake --profile slurm -s moments.smk --rerun-incomplete --retries 2