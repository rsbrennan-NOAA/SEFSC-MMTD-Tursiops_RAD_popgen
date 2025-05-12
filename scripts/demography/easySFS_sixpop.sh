#!/bin/bash
#SBATCH --job-name=easySFS_preview_sixpop
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --mem=20G
#SBATCH --partition=standard
#SBATCH --time=20:00:00

source ~/.bashrc
mamba activate easySFS

echo "six population start"

#remove hybs
# remove the hybrids from the clust file:
grep -v -f ~/Tursiops-RAD-popgen/analysis/pop_structure/newhybrids/hybrids.txt ~/Tursiops-RAD-popgen/analysis/pop_structure/sixpop.clust |\
     grep -v "Aduncus"  > ~/Tursiops-RAD-popgen/analysis/pop_structure/sixpop_subset_noHybrids.clust

~/bin/easySFS/easySFS.py  -i ~/Tursiops-RAD-popgen/analysis/variants/filtered.sixpop_ids.noMAF.vcf.gz -p ~/Tursiops-RAD-popgen/analysis/pop_structure/sixpop_subset_noHybrids.clust --preview -a

echo "six population done"

