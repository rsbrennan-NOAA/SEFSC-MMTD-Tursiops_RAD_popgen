#!/bin/bash
#SBATCH --job-name=SFS_c_six
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --mem=500G
#SBATCH --partition=himem
#SBATCH --time=7-00

source ~/.bashrc
mamba activate easySFS

cd ~/Tursiops-RAD-popgen/analysis/moments

echo "six population start"

~/bin/easySFS/easySFS.py  -i ~/Tursiops-RAD-popgen/analysis/variants/filtered.final.noMAF.vcf.gz -p ~/Tursiops-RAD-popgen/analysis/pop_structure/sixpop_all.clust -a -f -o ~/Tursiops-RAD-popgen/analysis/moments/sixpop_sfs --order Coastal_Atlantic,Coastal_Gulf,Intermediate_Atlantic,Intermediate_Gulf,Offshore_Atlantic,Offshore_Gulf --proj 224,64,58,20,40,50

echo "six population done"

