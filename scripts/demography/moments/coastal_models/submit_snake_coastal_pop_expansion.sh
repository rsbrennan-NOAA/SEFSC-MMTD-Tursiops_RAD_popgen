#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --partition standard

source ~/.bashrc
conda activate snakemake-9.6.2

cd /home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/coastal_models

for run_num in {7..12}; do
    echo "Starting run${run_num}"
    
    # Create directories
        mkdir -p ~/Tursiops-RAD-popgen/scripts/demography/moments//coastal_models/coastal_pop_expansion
        mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/coastal_pop_expansion/
    mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/coastal_pop_expansion/output_yaml
    
    # Unlock and run snakemake
    snakemake --profile slurm -s moments_coastal_pop_expansion.smk --unlock
    snakemake --profile slurm -s moments_coastal_pop_expansion.smk --rerun-incomplete --retries 2 --jobs 10 --config run_name="run${run_num}"
    
    echo "Completed run${run_num}"
done