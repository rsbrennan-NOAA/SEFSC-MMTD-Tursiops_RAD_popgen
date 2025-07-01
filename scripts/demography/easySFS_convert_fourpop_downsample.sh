#!/bin/bash
#SBATCH --job-name=SFS_c_four
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --mem=20G
#SBATCH --partition=standard
#SBATCH --time=7-00

source ~/.bashrc
mamba activate easySFS

cd ~/Tursiops-RAD-popgen/analysis/moments

grep -v -f ~/Tursiops-RAD-popgen/analysis/pop_structure/newhybrids/hybrids.txt ~/Tursiops-RAD-popgen/analysis/pop_structure/sixpop.clust |\
     grep -v "Aduncus"  > ~/Tursiops-RAD-popgen/analysis/pop_structure/sixpop_subset_noHybrids.clust

echo "four population start"

~/bin/easySFS/easySFS.py -i ~/Tursiops-RAD-popgen/analysis/variants/filtered.final.noMAF.vcf.gz -p ~/Tursiops-RAD-popgen/analysis/pop_structure/fourpop_all_noHybrids.clust -a -y -f -o ~/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs_downsample --order Coastal_Atlantic,Coastal_Gulf,Intermediate,Offshore -y --proj 40,40,40,40

echo "four population done"