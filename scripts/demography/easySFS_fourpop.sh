#!/bin/bash
#SBATCH --job-name=easySFS_preview_fourpop
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

# remove the hybrids from the clust file:
grep -v -f ~/Tursiops-RAD-popgen/analysis/pop_structure/newhybrids/hybrids.txt ~/Tursiops-RAD-popgen/analysis/pop_structure/fourpop.clust |\
     grep -v "Aduncus"  > ~/Tursiops-RAD-popgen/analysis/pop_structure/fourpop_subset_noHybrids.clust

echo "four population start"

~/bin/easySFS/easySFS.py  -i ~/Tursiops-RAD-popgen/analysis/variants/filtered.fourpop_ids.noMAF.vcf.gz -p ~/Tursiops-RAD-popgen/analysis/pop_structure/fourpop_subset_noHybrids.clust --preview -a

echo "four population done"

