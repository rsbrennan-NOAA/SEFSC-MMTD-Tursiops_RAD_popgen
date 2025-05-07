#!/bin/bash
#SBATCH --job-name=SFS_c_four
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --mem=210G
#SBATCH --partition=himem
#SBATCH --time=7-00

source ~/.bashrc
mamba activate easySFS

cd ~/Tursiops-RAD-popgen/analysis/moments

echo "four population start"

~/bin/easySFS/easySFS.py  -i ~/Tursiops-RAD-popgen/analysis/variants/filtered.final.noMAF.vcf.gz -p ~/Tursiops-RAD-popgen/analysis/pop_structure/fourpop_all_noHybrids.clust -a -f -o ~/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs --order Coastal_Atlantic,Coastal_Gulf,Intermediate,Offshore  --proj 224,64,66,66

echo "four population done"
