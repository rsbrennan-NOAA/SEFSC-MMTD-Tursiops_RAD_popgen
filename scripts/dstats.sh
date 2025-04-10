#!/bin/bash


cd ~/Tursiops-RAD-popgen/analysis/pop_structure/dstats

#Aduncus	SRR5357657
#Aduncus	SRR5357656
#Aduncus	SRR5357655

~/bin/Dsuite/Build/Dsuite Dtrios --out-prefix fourpop_all --no-f4-ratio ~/Tursiops-RAD-popgen/analysis/variants/tursiops_aduncus_LDthin.vcf.gz fourpop_all_dstats.pop

~/bin/Dsuite/Build/Dsuite Dtrios --out-prefix sixpop_all --no-f4-ratio ~/Tursiops-RAD-popgen/analysis/variants/tursiops_aduncus_LDthin.vcf.gz sixpop_all_dstats.pop

cut -f 2 fourpop_all_dstats.pop | sort|uniq | grep -v Outgroup > plot_order_fourpop.txt
cut -f 2 sixpop.pop | sort| uniq |grep -v Outgroup > plot_order_sixpop.txt
