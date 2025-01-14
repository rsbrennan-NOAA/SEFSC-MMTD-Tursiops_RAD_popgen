#!/bin/bash
#SBATCH --job-name=filter2
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=40:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16

INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants

vcftools --gzvcf ${INDIR}/filtered.6.vcf.gz --max-meanDP 124 --exclude-positions ~/Tursiops-NC-PopulationAssignment-RAD/scripts/HD_exclude.txt --recode --recode-INFO-all --stdout |  bgzip > ${INDIR}/filtered.final.vcf.gz
