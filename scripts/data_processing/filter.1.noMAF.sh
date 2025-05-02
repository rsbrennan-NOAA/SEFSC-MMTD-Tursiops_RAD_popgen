#!/bin/bash
#SBATCH --job-name=filter1
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=16:00:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16
module load bio/bcftools/1.11
module load bio/htslib/1.19

INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants
OUTDIR=~/Tursiops-RAD-popgen/analysis/variants

# note, starting with existing files bc nothing changes before this step.

echo "start 2nd filtering"
vcftools --gzvcf ${INDIR}/filtered.2.vcf.gz \
	        --minDP 10 --remove-indels --max-alleles 2 --min-alleles 2 --minQ 30  \
		--recode --recode-INFO-all --stdout | \
		bgzip > ${OUTDIR}/filtered.3.noMAF.vcf.gz

echo "finish second filtering"

NUM_VARIANTS3=$(zcat ${OUTDIR}/filtered.3.noMAF.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after second filters: ${NUM_VARIANTS3}"

# check missingness:
cd ${INDIR}
vcftools --gzvcf ${OUTDIR}/filtered.3.noMAF.vcf.gz --missing-indv


