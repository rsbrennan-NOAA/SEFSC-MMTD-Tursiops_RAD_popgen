#!/bin/bash
module load bio/plink/1.90b6.23

plink --vcf ~/Tursiops-RAD-popgen/analysis/variants/LDthin_numCorrect.vcf \
	--distance square --out ~/Tursiops-RAD-popgen/analysis/pop_structure/LDthin \
	--allow-extra-chr --double-id
