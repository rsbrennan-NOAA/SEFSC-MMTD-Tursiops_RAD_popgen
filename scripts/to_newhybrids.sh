#!/bin/bash
#SBATCH --job-name=newhybrids
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --time=1:00:00


module load bio/vcftools/0.1.16
module load bio/htslib/1.19
module load bio/plink/1.90b6.23

cd ~/Tursiops-RAD-popgen/analysis/variants

vcftools --gzvcf ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_ids_LDthin.vcf.gz --snps ~/Tursiops-RAD-popgen/analysis/snp_set_gulfVSoffshore.txt --recode --recode-INFO-all --stdout | bgzip > ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_ids_LDthin_gulfVSoffshore.vcf.gz

plink --vcf ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_ids_LDthin_gulfVSoffshore.vcf.gz --recode --allele1234 --allow-extra-chr --out all_variants
#65 variants and 345 people pass filters and QC.

cut -f 7- -d " " all_variants.ped > all_variants.1.ped
sed -i 's/0/-1/g' all_variants.1.ped
awk '{print NR, $0}' all_variants.1.ped > all_variants.2.ped
#head all_variants.2.ped | cut -f 1-10 -d " "
# Count number of individuals in the PED file
num_individuals=$(wc -l < all_variants.1.ped)
# Count number of loci in the MAP file
num_loci=$(wc -l < all_variants.map)
# get loci
loci=$(awk '{printf "%s ", $2}' all_variants.map | sed 's/ $//')
echo $loci
# Create a header file with the extracted values
echo "NumIndivs $num_individuals" > temp_header
echo "NumLoci $num_loci" >> temp_header
echo "Digits 1" >> temp_header
echo "Format NonLumped" >> temp_header
# Create the locus names line
echo "LocusNames $(echo $loci | tr '\n' ' ')" >> temp_header
cat all_variants.2.ped >> temp_header
mv temp_header all_variants.3.ped
#head ~/Tursiops-RAD-popgen/analysis/variants/all_variants.3.ped  | cut -f 1-10 -d " "


