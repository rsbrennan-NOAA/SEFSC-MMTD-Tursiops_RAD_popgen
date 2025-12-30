#!/bin/bash
#SBATCH --job-name=dstats
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --time=1:00:00

cd ~/Tursiops-RAD-popgen/analysis/pop_structure/dstats


#Aduncus	SRR5357657
#Aduncus	SRR5357656
#Aduncus	SRR5357655

~/bin/Dsuite/Build/Dsuite Dtrios --out-prefix fourpop_all --no-f4-ratio ~/Tursiops-RAD-popgen/analysis/variants/tursiops_aduncus_LDthin.vcf.gz fourpop_all_dstats.pop

~/bin/Dsuite/Build/Dsuite Dtrios --out-prefix sixpop_all --no-f4-ratio ~/Tursiops-RAD-popgen/analysis/variants/tursiops_aduncus_LDthin.vcf.gz sixpop_all_dstats.pop

cut -f 2 fourpop_all_dstats.pop | sort|uniq | grep -v Outgroup > plot_order_fourpop.txt
cut -f 2 sixpop.pop | sort| uniq |grep -v Outgroup > plot_order_sixpop.txt

# repeat with no hybrids
~/bin/Dsuite/Build/Dsuite Dtrios --out-prefix fourpop_nohybs --no-f4-ratio ~/Tursiops-RAD-popgen/analysis/variants/tursiops_aduncus_LDthin.vcf.gz fourpop_nohybs_dstats.pop


