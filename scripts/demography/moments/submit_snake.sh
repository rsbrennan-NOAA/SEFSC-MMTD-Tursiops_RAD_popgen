#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --partition standard

source ~/.bashrc
conda activate snakemake-9.6.2

cd /home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments


# expansion run 4:

mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_expansion/
mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_expansion/output_yaml

snakemake --profile slurm -s moments.smk --unlock

snakemake --profile slurm -s moments.smk --rerun-incomplete --retries 2 --config run_name="run_4"

# expansion run 4:

mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_expansion/
mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_expansion/output_yaml

snakemake --profile slurm -s moments.smk --unlock

snakemake --profile slurm -s moments.smk --rerun-incomplete --retries 2 --config run_name="run_5"

# simple run 1

mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/
mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/output_yaml

snakemake --profile slurm -s moments_mod7_simple.smk --unlock

snakemake --profile slurm -s moments_mod7_simple.smk --rerun-incomplete --retries 2 --config run_name="run_1"


# simple run2

mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/
mkdir -p /home/rbrennan/Tursiops-RAD-popgen/analysis/moments/mod7_simple/output_yaml

snakemake --profile slurm -s moments_mod7_simple.smk --unlock

snakemake --profile slurm -s moments_mod7_simple.smk --rerun-incomplete --retries 2 --config run_name="run_2"
