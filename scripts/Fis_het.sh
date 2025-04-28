#!/bin/bash
#SBATCH --job-name=populations
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -c 14
#SBATCH --mem=20G
#SBATCH --partition=standard
#SBATCH --time=2:00:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16
module load bio/stacks/2.65

populations --in-vcf ~/Tursiops-RAD-popgen/analysis/variants/filtered.final.vcf.gz --popmap ~/Tursiops-RAD-popgen/analysis/pop_structure/fourpop.clust --fstats --threads 14 --out-path ~/Tursiops-RAD-popgen/analysis/diversity/populations_fourpops

populations --in-vcf ~/Tursiops-RAD-popgen/analysis/variants/filtered.final.vcf.gz --popmap ~/Tursiops-RAD-popgen/analysis/pop_structure/sixpop.clust --fstats --threads 14 --out-path ~/Tursiops-RAD-popgen/analysis/diversity/populations_sixpops

# remove hybrids

vcftools --gzvcf ~/Tursiops-RAD-popgen/analysis/variants/filtered.final.vcf.gz --not-chr NC_047055.1 --remove ~/Tursiops-RAD-popgen/analysis/pop_structure/newhybrids/hybrids.txt --recode --recode-INFO-all --stdout |  bgzip > ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_nohybrids.vcf.gz

populations --in-vcf ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_nohybrids.vcf.gz --popmap ~/Tursiops-RAD-popgen/analysis/pop_structure/fourpop.clust --fstats --threads 14 --out-path ~/Tursiops-RAD-popgen/analysis/diversity/populations_noHybs_fourpops

populations --in-vcf ~/Tursiops-RAD-popgen/analysis/variants/filtered.final_nohybrids.vcf.gz --popmap ~/Tursiops-RAD-popgen/analysis/pop_structure/sixpop.clust --fstats --threads 14 --out-path ~/Tursiops-RAD-popgen/analysis/diversity/populations_noHybs_sixpops

