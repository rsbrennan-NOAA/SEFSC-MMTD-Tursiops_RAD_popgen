#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -c 3
#SBATCH --mem=10G
#SBATCH --partition=standard
#SBATCH --time=6:00:00

source ~/.bashrc

# load fastqc and multiqc
mamba activate multiqc-1.17

module load bio/fastqc/0.11.9


cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus

mkdir -p fastqc_rawData
fastqc -t 14 --noextract -o fastqc_rawData/ *.fastq.gz

multiqc fastqc_rawData/ -o  fastqc_rawData/ --filename multiqc_report_rawData
