#!/bin/bash
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --mail-type=END
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --partition=standard
#SBATCH --mem=20G
#SBATCH --time=05:00:00
#SBATCH --job-name=bam_stats
#SBATCH --output=%x.%A.%a.out
#SBATCH --error=%x.%A.%a.err

module load bio/samtools
cd ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus/aligned/mkdup

echo "sample,total_reads,total_mapped,map_percent,duplicates,duplicate_percent,mapped_q20,map_percent_q20"> ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus/aligned/mkdup/count.aligned.txt

for i in $(ls *.bam | cut -f -1 -d "." | uniq );do

samtools flagstat ${i}.bam > tmp.${i}.txt

TOTAL=$(grep "in total" tmp.${i}.txt | cut -f 1 -d " ")
MAPPED=$(grep "mapped (" tmp.${i}.txt | cut -f 1 -d " " | head -n 1)
PERCENT=$(grep "mapped (" tmp.${i}.txt | cut -f 2 -d "(" | cut -f 1 -d "%" | head -n 1)
MAPPED_q=$(samtools view -F 4 -q 20 ${i}.bam | wc -l)
PERCENT_q=$(echo "scale=2 ; $MAPPED_q / $TOTAL" | bc)
DUPS=$(grep "duplicates" tmp.${i}.txt | cut -f 1 -d " ")
DUP_PERCENT=$(echo "scale=2 ; $DUPS  / $TOTAL" | bc)

echo "${i},$TOTAL,$MAPPED,$PERCENT,$DUPS,$DUP_PERCENT,$MAPPED_q,$PERCENT_q"

done  >> ~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus/aligned/mkdup/count.aligned.txt
