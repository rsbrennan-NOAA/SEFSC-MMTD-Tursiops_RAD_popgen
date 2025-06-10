#!/bin/bash
#SBATCH --job-name=vcf_concat
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --partition=standard
#SBATCH --time=72:00:00

# load modules
module load bio/bcftools/1.11
module load bio/htslib/1.19
OUTDIR=~/Tursiops-RAD-popgen/analysis/variants/aduncus/
INDIR=~/Tursiops-RAD-popgen/analysis/variants/aduncus/chromosomes

# make list of VCF files
ls $INDIR/variants_raw_*.vcf.gz > $OUTDIR/vcf_list.txt

echo "start merge"

# merge vcfs
#bcftools concat -f $OUTDIR/vcf_list.txt -O z -o $OUTDIR/variants_raw_merged.vcf.gz

# index
#tabix -p vcf $OUTDIR/variants_raw_merged.vcf.gz

echo "done with merge"

#
# filter -----------------------------------------

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16
module load bio/bcftools/1.11
module load bio/htslib/1.19

INDIR=~/Tursiops-RAD-popgen/analysis/variants/aduncus/

echo "starting first missingness filter"

bcftools view -i 'F_MISSING <= 0.4 && FMT/DP >= 5' -Oz -o ${INDIR}/filtered.1.vcf.gz ${INDIR}/variants_raw_merged.vcf.gz

# Index the output
bcftools index -t ${INDIR}/filtered.1.vcf.gz

echo "done with first missingness filter"

NUM_VARIANTS1=$(zcat ${INDIR}/filtered.1.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after first missingness filter: ${NUM_VARIANTS1}"

echo "running vcfallelicprimitives"
vcfallelicprimitives ${INDIR}/filtered.1.vcf.gz --keep-info --keep-geno | \
	        vcfstreamsort |  bgzip > ${INDIR}/filtered.2.vcf.gz

echo "done with allelic primitive"
NUM_VARIANTS2=$(zcat ${INDIR}/filtered.2.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after primitives step: ${NUM_VARIANTS2}"

echo "start 2nd filtering"

vcftools --gzvcf ${INDIR}/filtered.2.vcf.gz \
    --maxDP 800 \
    --remove-indels \
    --max-alleles 2 \
    --min-alleles 2 \
    --minQ 20 \
    --recode --recode-INFO-all --stdout | \
    bgzip > ${INDIR}/aduncus.vcf.gz

echo "finish second filtering"

NUM_VARIANTS3=$(zcat ${INDIR}/aduncus.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after second filters: ${NUM_VARIANTS3}"

bcftools index -t ${INDIR}/aduncus.vcf.gz 
