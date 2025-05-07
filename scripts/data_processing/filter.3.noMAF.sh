#!/bin/bash
#SBATCH --job-name=filter3
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

vcftools --gzvcf ${INDIR}/filtered.6.noMAF.vcf.gz --max-meanDP 136 --exclude-positions ~/Tursiops-RAD-popgen/scripts/HD_exclude_noMAF.txt --recode --recode-INFO-all --stdout |  bgzip > ${INDIR}/filtered.final_MitoIncluded.noMAF.vcf.gz

# remove the mitochondrial chromosome

vcftools --gzvcf ${INDIR}/filtered.final_MitoIncluded.noMAF.vcf.gz \
        --chr NC_047034.1 \
        --chr NC_047035.1 \
        --chr NC_047036.1 \
        --chr NC_047037.1 \
        --chr NC_047038.1 \
        --chr NC_047039.1 \
        --chr NC_047040.1 \
        --chr NC_047041.1 \
        --chr NC_047042.1 \
        --chr NC_047043.1 \
        --chr NC_047044.1 \
        --chr NC_047045.1 \
        --chr NC_047046.1 \
        --chr NC_047047.1 \
        --chr NC_047048.1 \
        --chr NC_047049.1 \
        --chr NC_047050.1 \
        --chr NC_047051.1 \
        --chr NC_047052.1 \
        --chr NC_047053.1 \
        --chr NC_047054.1 \
        --chr NC_047055.1 \
	--remove ~/Tursiops-RAD-popgen/analysis/pop_structure/newhybrids/hybrids.txt \
	--recode --recode-INFO-all --stdout |  bgzip > ${INDIR}/filtered.final.noMAF.vcf.gz

bcftools annotate --set-id '%CHROM\_%POS'  ${INDIR}/filtered.final.noMAF.vcf.gz  -Oz -o ${INDIR}/filtered.final_ids.noMAF.vcf.gz

bcftools index -t ${INDIR}/filtered.final_MitoIncluded.noMAF.vcf.gz
bcftools index -t ${INDIR}/filtered.final.noMAF.vcf.gz
bcftools index -t ${INDIR}/filtered.final_ids.noMAF.vcf.gz

