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
#bcftools query -f '%CHROM\t%POS\n' filtered_lenient.final_ids_LDthin.vcf.gz > rad_positions.txt

# keep only rad positions
#bcftools view -R rad_positions.txt ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus.vcf.gz |bgzip > ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz
#tabix ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz

zcat ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz | grep -v "^#" | wc -l
# 2364

#merge

#bcftools merge ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz filtered_lenient.final_ids_LDthin.vcf.gz| bgzip > tursiops_aduncus_1.vcf.gz
#tabix tursiops_aduncus_1.vcf.gz

#zcat tursiops_aduncus_1.vcf.gz | grep -v '^#' | wc -l
# 7557

# keep only the ones with sites in aduncus:
#bcftools query -f '%CHROM\t%POS\n' ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz > rad_positions_aduncus.txt
#bcftools view -R rad_positions_aduncus.txt tursiops_aduncus_1.vcf.gz |bgzip > tursiops_aduncus_lenient_noMissing.vcf.gz
#zcat tursiops_aduncus_noMissing.vcf.gz | grep -v '^#' | wc -l
#cp tursiops_aduncus_1.vcf.gz tursiops_aduncus_lenient_All.vcf.gz 
# 

#bcftools index -t tursiops_aduncus_LDthin.vcf.gz


######
# original set:


# get rad positions
bcftools query -f '%CHROM\t%POS\n' filtered.final_ids.vcf.gz > rad_positions.txt
cat rad_positions.txt | wc -l
# 4495

# keep only rad positions
bcftools view -R rad_positions.txt ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus.vcf.gz |bgzip > ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz
tabix ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz

zcat ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz | grep -v "^#" | wc -l
# 2364

#merge

bcftools merge ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz filtered.final_ids.vcf.gz| bgzip > tursiops_aduncus_1.vcf.gz
tabix tursiops_aduncus_1.vcf.gz

zcat tursiops_aduncus_All.vcf.gz | grep -v '^#' | wc -l
# 4495

# keep only the ones with sites in aduncus:
bcftools query -f '%CHROM\t%POS\n' ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites.vcf.gz > rad_positions_aduncus.txt
bcftools view -R rad_positions_aduncus.txt tursiops_aduncus_1.vcf.gz |bgzip > tursiops_aduncus.vcf.gz
zcat tursiops_aduncus.vcf.gz | grep -v '^#' | wc -l

bcftools index -t tursiops_aduncus.vcf.gz

plink --vcf tursiops_aduncus.vcf.gz \
        --indep-pairwise 50 5 0.2 --allow-extra-chr --double-id \
        --out variants_tursiops_aduncus_pruned
#Pruning complete.  1505 of 2795 variants removed.
# make thinned vcf
vcftools --gzvcf tursiops_aduncus.vcf.gz --snps variants_tursiops_aduncus_pruned.prune.in  --recode --recode-INFO-all --stdout | \
        bgzip > tursiops_aduncus_LDthin.vcf.gz

bcftools index -t tursiops_aduncus_LDthin.vcf.gz


#######################################################
# remove x

bcftools query -f '%CHROM\t%POS\n' filtered.final_ids_noX.vcf.gz > rad_positions_noX.txt
bcftools view -R rad_positions_noX.txt ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus.vcf.gz |bgzip > ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites_noX.vcf.gz
tabix ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites_noX.vcf.gz

bcftools merge ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites_noX.vcf.gz filtered.final_ids_noX.vcf.gz | bgzip > tursiops_aduncus_noX_1.vcf.gz
tabix tursiops_aduncus_noX_1.vcf.gz
zcat tursiops_aduncus_1.vcf.gz | grep -v '^#' | wc -l

# keep only the ones with sites in aduncus:
bcftools query -f '%CHROM\t%POS\n' ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites_noX.vcf.gz > rad_positions_aduncus_noX.txt
bcftools view -R rad_positions_aduncus_noX.txt tursiops_aduncus_noX_1.vcf.gz |bgzip > tursiops_aduncus_noX.vcf.gz
tabix tursiops_aduncus_noX.vcf.gz
zcat tursiops_aduncus_noX.vcf.gz | grep -v '^#' | wc -l


plink --vcf tursiops_aduncus_noX.vcf.gz \
        --indep-pairwise 50 5 0.2 --allow-extra-chr --double-id \
        --out variants_tursiops_aduncus_noX_pruned
#Pruning complete. 
# make thinned vcf
vcftools --gzvcf tursiops_aduncus_noX.vcf.gz --snps variants_tursiops_aduncus_noX_pruned.prune.in  --recode --recode-INFO-all --stdout | \
        bgzip > tursiops_aduncus_noX_LDthin.vcf.gz

bcftools index -t tursiops_aduncus_noX_LDthin.vcf.gz




############---------------------------------------------------
# all above are ld thin loci. Do a set with all loci:

#bcftools query -f '%CHROM\t%POS\n' filtered.final.vcf.gz > rad_positions_nothin.txt

#bcftools view -R rad_positions_nothin.txt ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus.vcf.gz |bgzip > ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites_nothin.vcf.gz


#tabix ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites_nothin.vcf.gz

#zcat ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites_nothin.vcf.gz | grep -v "^#" | wc -l
#2901 

#merge

#bcftools merge ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites_nothin.vcf.gz filtered.final.vcf.gz| bgzip > tursiops_aduncus_nothin.1.vcf.gz
#tabix tursiops_aduncus_nothin.1.vcf.gz

#zcat tursiops_aduncus_nothin.1.vcf.gz | grep -v '^#' | wc -l
# 7818

# keep only the ones with sites in aduncus:
#bcftools query -f '%CHROM\t%POS\n' ~/Tursiops-RAD-popgen/analysis/variants/aduncus/aduncus_radsites_nothin.vcf.gz > rad_positions_aduncus_nothin.txt
#bcftools view -R rad_positions_aduncus_nothin.txt tursiops_aduncus_nothin.1.vcf.gz |bgzip > tursiops_aduncus_nothin.vcf.gz
#zcat tursiops_aduncus_nothin.vcf.gz | grep -v '^#' | wc -l
# 2901 

#bcftools index -t tursiops_aduncus_nothin.vcf.gz










