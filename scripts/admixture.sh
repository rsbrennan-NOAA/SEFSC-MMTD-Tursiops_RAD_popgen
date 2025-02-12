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


module load  bio/admixture/1.3.0
module load bio/vcftools/0.1.16
module load bio/plink/1.90b6.23

cd ~/Tursiops-RAD-popgen/analysis/pop_structure

zcat ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_ids_LDthin.vcf.gz| \
grep -v "^##contig=<ID=N[CW]_"  | awk -F'\t' 'BEGIN {OFS="\t"} 
    /^#/ {print; next} 
    {split($1, a, /\./); sub(/^N[CW]_/, "", a[1]); $1 = a[1]; print}' > ~/Tursiops-RAD-popgen/analysis/variants/LDthin_numCorrect.vcf

plink --vcf ~/Tursiops-RAD-popgen/analysis/variants/LDthin_numCorrect.vcf \
--make-bed --out ~/Tursiops-RAD-popgen/analysis/variants/LDthin_numCorrect \
--allow-extra-chr --double-id


for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; \
do admixture -j4 -s time --cv ~/Tursiops-RAD-popgen/analysis/variants/LDthin_numCorrect.bed $K | tee log${K}.out; done

grep -h CV log*.out | cut -f 3- -d " " > cv.txt




