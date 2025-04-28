#!/bin/bash
#SBATCH --job-name=ld_thin
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=01:00:00


module load bio/vcftools/0.1.16
module load bio/plink/1.90b6.23
module load bio/htslib/1.19
module load bio/bcftools/1.11

INDIR=~/Tursiops-RAD-popgen/analysis/variants

plink --vcf ${INDIR}/filtered.final_ids.vcf.gz \
	        --indep-pairwise 50 5 0.2 --allow-extra-chr --double-id \
		        --out variants_pruned

# make thinned vcf
vcftools --gzvcf ${INDIR}/filtered.final_ids.vcf.gz --snps variants_pruned.prune.in  --recode --recode-INFO-all --stdout | \
	        bgzip > ${INDIR}/filtered.final_ids_LDthin.vcf.gz

bcftools index -t ${INDIR}/filtered.final_ids_LDthin.vcf.gz
