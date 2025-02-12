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

# subset with subset.txt

plink --vcf ~/Tursiops-RAD-popgen/analysis/variants/LDthin_numCorrect.vcf \
	--keep subset.txt \
	--make-bed --out ~/Tursiops-RAD-popgen/analysis/variants/LDthin_numCorrect_subset \
	--allow-extra-chr --double-id


for K in 1 2 3 4 5 6 7 8 9 10; \
do admixture -j4 -s time --cv ~/Tursiops-RAD-popgen/analysis/variants/LDthin_numCorrect_subset.bed $K | tee log_subset.${K}.out; done

grep -h CV log_subset.*.out | cut -f 3- -d " " > cv_subset.txt


