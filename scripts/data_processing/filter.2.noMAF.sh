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
#SBATCH --time=3:00:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16

INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants
OUTDIR=~/Tursiops-RAD-popgen/analysis/variants/

vcftools --gzvcf ${INDIR}/filtered.3.vcf.gz --remove ~/Tursiops-RAD-popgen/scripts/rm_missing.txt  --recode --recode-INFO-all --stdout |  bgzip > ${OUTDIR}/filtered.4.noMAF.vcf.gz 

vcftools --gzvcf ${OUTDIR}/filtered.4.vcf.gz  --max-missing 0.7 --min-meanDP 10 --recode --recode-INFO-all  --stdout  >  ${OUTDIR}/filtered.5.noMAF.vcf

NUM_VARIANTS1=$(cat ${OUTDIR}/filtered.5.noMAF.vcf | grep -v '^#' | wc -l)
echo "Number of variants after missing, no maf filters: ${NUM_VARIANTS1}"

vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" ${OUTDIR}/filtered.5.noMAF.vcf |bgzip > ${OUTDIR}/filtered.6.noMAF.vcf.gz

NUM_VARIANTS2=$(zcat ${OUTDIR}/filtered.6.noMAF.vcf.gz | grep -v '^#' | wc -l)
 echo "Number of variants after AB  filters: ${NUM_VARIANTS2}"

#depths
vcftools --gzvcf ${OUTDIR}/filtered.6.noMAF.vcf.gz --site-mean-depth --out ${OUTDIR}/filtered.6.depth
