#!/bin/bash
#SBATCH --job-name=easySFS_preview_sixpop
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --mem=380G
#SBATCH --partition=himem
#SBATCH --time=20:00:00

source ~/.bashrc
mamba activate easySFS

echo "six population start"

~/bin/easySFS/easySFS.py  -i ~/Tursiops-RAD-popgen/analysis/variants/filtered.final.noMAF.vcf.gz -p ~/Tursiops-RAD-popgen/analysis/pop_structure/sixpop_all.clust --preview -a

echo "six population done"

