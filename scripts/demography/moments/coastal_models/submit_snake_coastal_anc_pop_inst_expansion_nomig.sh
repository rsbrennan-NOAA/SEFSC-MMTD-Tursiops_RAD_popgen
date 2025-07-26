#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --partition standard

source ~/.bashrc
conda activate snakemake-9.6.2

MODEL="coastal_anc_pop_inst_expansion_nomig"

cd /home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/coastal_models
for run_num in {1..20}; do
    echo "Starting run${run_num}"
    
    # Create directories
        mkdir -p ~/Tursiops-RAD-popgen/scripts/demography/moments/coastal_models/${MODEL}
        mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/${MODEL}/
    mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/${MODEL}/output_yaml
    
    # Unlock and run snakemake
    snakemake --profile slurm -s moments_coastal_generic.smk --unlock
    snakemake --profile slurm -s moments_coastal_generic.smk --rerun-incomplete --retries 2 --jobs 10 \
        --config model_name="${MODEL}" run_name="run${run_num}"
    
    echo "Completed ${MODEL} run${run_num}"
done