#!/bin/bash
#SBATCH --job-name=filter2
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=40:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16
module load bio/bcftools/1.11

INDIR=~/Tursiops-RAD-popgen/analysis/variants

vcftools --gzvcf ${INDIR}/filtered.6.vcf.gz --max-meanDP 136 --exclude-positions ~/Tursiops-RAD-popgen/scripts/HD_exclude.txt --recode --recode-INFO-all --stdout |  bgzip > ${INDIR}/filtered.final_MitoIncluded.vcf.gz

# remove the mitochondrial chromosome

vcftools --gzvcf ${INDIR}/filtered.final_MitoIncluded.vcf.gz --not-chr NC_012059.1 --not-chr NW_022983433.1 --recode --recode-INFO-all --stdout |  bgzip > ${INDIR}/filtered.final.vcf.gz

bcftools annotate --set-id '%CHROM\_%POS'  ${INDIR}/filtered.final.vcf.gz  -Oz -o ${INDIR}/filtered.final_ids.vcf.gz

bcftools index -t ${INDIR}/filtered.final_MitoIncluded.vcf.gz
bcftools index -t ${INDIR}/filtered.final.vcf.gz
bcftools index -t ${INDIR}/filtered.final_ids.vcf.gz




