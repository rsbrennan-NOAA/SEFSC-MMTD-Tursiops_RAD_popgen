#!/bin/bash
#SBATCH --job-name=intersect
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --time=1:00:00

module load bio/vcftools/0.1.16
module load bio/plink/1.90b6.23
module load bio/htslib/1.19
module load bio/bcftools/1.11

cd ~/Tursiops-RAD-popgen/analysis/variants

# get rad positions
bcftools query -f '%CHROM\t%POS\n' filtered.final_ids.vcf.gz > rad_positions.txt
# keep only rad positions
bcftools view -R rad_positions.txt ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus.vcf.gz |bgzip > ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz
tabix ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz
#merge
bcftools merge ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz filtered.final_ids.vcf.gz|bgzip > tursiops_aduncus_1.vcf.gz
tabix tursiops_aduncus_1.vcf.gz

zcat tursiops_aduncus_1.vcf.gz | grep -v '^#' | wc -l

# keep only the ones with sites in aduncus:
bcftools query -f '%CHROM\t%POS\n' ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz > rad_positions_aduncus.txt
bcftools view -R rad_positions_aduncus.txt tursiops_aduncus_1.vcf.gz |bgzip > tursiops_aduncus.vcf.gz
tabix tursiops_aduncus.vcf.gz
zcat tursiops_aduncus.vcf.gz | grep -v '^#' | wc -l
#2795

#LD thin:

plink --vcf tursiops_aduncus.vcf.gz \
        --indep-pairwise 50 5 0.2 --allow-extra-chr --double-id \
        --out variants_tursiops_aduncus_pruned
#Pruning complete.  1505 of 2795 variants removed.
# make thinned vcf
vcftools --gzvcf tursiops_aduncus.vcf.gz --snps variants_tursiops_aduncus_pruned.prune.in  --recode --recode-INFO-all --stdout | \
        bgzip > tursiops_aduncus_LDthin.vcf.gz

bcftools index -t tursiops_aduncus_LDthin.vcf.gz
