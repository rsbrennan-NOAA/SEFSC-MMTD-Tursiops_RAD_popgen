#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -c 4
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=03:00:00


source ~/.bashrc

mamba activate vcflib-1.0.9
module load  bio/admixture/1.3.0
module load bio/vcftools/0.1.16
module load bio/plink/1.90b6.23

cd ~/Tursiops-RAD-popgen/analysis/pop_structure

# drop sex chr
vcftools --gzvcf ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_ids_LDthin.vcf.gz --not-chr NC_047055.1 --recode --recode-INFO-all --stdout |  bgzip > ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_ids_LDthin_noX.vcf.gz


zcat ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_ids_LDthin_noX.vcf.gz| \
grep -v "^##contig=<ID=N[CW]_"  | awk -F'\t' 'BEGIN {OFS="\t"} 
    /^#/ {print; next} 
    {split($1, a, /\./); sub(/^N[CW]_/, "", a[1]); $1 = a[1]; print}' > ~/Tursiops-RAD-popgen/analysis/variants/LDthin_noX_numCorrect.vcf

plink --vcf ~/Tursiops-RAD-popgen/analysis/variants/LDthin_noX_numCorrect.vcf \
--make-bed --out ~/Tursiops-RAD-popgen/analysis/variants/LDthin_noX_numCorrect \
--allow-extra-chr --double-id


rm log*.out

for K in 1 2 3 4 5 6 7 8 9 10; \
do admixture -j4 -s time --cv ~/Tursiops-RAD-popgen/analysis/variants/LDthin_noX_numCorrect.bed $K | tee log${K}.out; done


grep -h CV log*.out | cut -f 3- -d " " > cv.txt




