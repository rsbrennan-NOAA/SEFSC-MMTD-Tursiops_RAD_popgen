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

source ~/.bashrc

mamba activate vcflib-1.0.9

module load bio/vcftools/0.1.16
module load bio/plink/1.90b6.23
module load bio/htslib/1.19
module load bio/bcftools/1.11

cd ~/Tursiops-RAD-popgen/analysis/variants

# remove x chrom

#vcftools --gzvcf ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_ids.vcf.gz --not-chr NC_047055.1 --recode --recode-INFO-all --stdout |  bgzip > ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_ids_noX.vcf.gz

#tabix filtered.final_ids_noX.vcf.gz

# get rad positions
bcftools query -f '%CHROM\t%POS\n' filtered_lenient.final_ids_LDthin.vcf.gz > rad_positions.txt

# keep only rad positions
bcftools view -R rad_positions.txt ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus.vcf.gz |bgzip > ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz
tabix ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz

zcat ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz | grep -v "^#" | wc -l
# 2364

#merge

bcftools merge ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz filtered_lenient.final_ids_LDthin.vcf.gz| bgzip > tursiops_aduncus_1.vcf.gz
tabix tursiops_aduncus_1.vcf.gz

zcat tursiops_aduncus_1.vcf.gz | grep -v '^#' | wc -l
# 7557

# keep only the ones with sites in aduncus:
bcftools query -f '%CHROM\t%POS\n' ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz > rad_positions_aduncus.txt
bcftools view -R rad_positions_aduncus.txt tursiops_aduncus_1.vcf.gz |bgzip > tursiops_aduncus_lenient_noMissing.vcf.gz
zcat tursiops_aduncus_noMissing.vcf.gz | grep -v '^#' | wc -l
cp tursiops_aduncus_1.vcf.gz tursiops_aduncus_lenient_All.vcf.gz 
# 

bcftools index -t tursiops_aduncus_LDthin.vcf.gz


######
# original set:


# get rad positions
bcftools query -f '%CHROM\t%POS\n' filtered.final_ids_LDthin.vcf.gz > rad_positions.txt
cat rad_positions.txt | wc -l
# 4495

# keep only rad positions
bcftools view -R rad_positions.txt ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus.vcf.gz |bgzip > ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz
tabix ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz

zcat ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz | grep -v "^#" | wc -l
# 2364

#merge

bcftools merge ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz filtered.final_ids_LDthin.vcf.gz| bgzip > tursiops_aduncus_All.vcf.gz
tabix tursiops_aduncus_All.vcf.gz

zcat tursiops_aduncus_All.vcf.gz | grep -v '^#' | wc -l
# 4495

# keep only the ones with sites in aduncus:
bcftools query -f '%CHROM\t%POS\n' ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz > rad_positions_aduncus.txt
bcftools view -R rad_positions_aduncus.txt tursiops_aduncus_All.vcf.gz |bgzip > tursiops_aduncus_noMissing.vcf.gz
zcat tursiops_aduncus_noMissing.vcf.gz | grep -v '^#' | wc -l
cp tursiops_aduncus_1.vcf.gz tursiops_aduncus_lenient_All.vcf.gz 
# 

bcftools index -t tursiops_aduncus_LDthin.vcf.gz

