#!/bin/bash
#SBATCH --job-name=filter3
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH --partition=standard
#SBATCH --time=8:00:00

source ~/.bashrc

module load bio/sratoolkit/3.0.7 
mamba activate vcflib-1.0.9

cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus

prefetch SRX2653497 --max-size 60GB
fasterq-dump SRR5357657  --split-files --outdir ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus
bgzip --threads 8 SRR5357657_1.fastq
bgzip --threads 8 SRR5357657_2.fastq

prefetch SRX2653496 --max-size 90GB
fasterq-dump SRR5357656 --split-files --outdir ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus
bgzip --threads 8 SRR5357656_1.fastq
bgzip --threads 8 SRR5357656_2.fastq

prefetch SRX2653495 --max-size 90GB
fasterq-dump SRR5357655 --split-files --outdir ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus
bgzip --threads 8 SRR5357655_1.fastq
bgzip --threads 8 SRR5357655_2.fastq
