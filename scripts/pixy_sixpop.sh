#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 14
#SBATCH --mem=60G
#SBATCH --partition=standard
#SBATCH --time=20:00:00

source ~/.bashrc

mamba activate pixy

cd ~/Tursiops-RAD-popgen/analysis/diversity/
POP=sixpop

pixy --stats pi \
--vcf ../variants/combined_filtered_invariant.vcf.gz \
--populations ~/Tursiops-RAD-popgen/analysis/pop_structure/${POP}.clust \
--window_size 1000 \
--n_cores 14 \
--output_prefix 1kb_${POP}
