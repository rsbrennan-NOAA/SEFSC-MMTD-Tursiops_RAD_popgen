
module load bio/vcftools/0.1.16
module load bio/plink/1.90b6.23
module load bio/htslib/1.19
module load bio/bcftools/1.11

cd ~/Tursiops-RAD-popgen/analysis/variants
vcftools --gzvcf $filtered.4.vcf.gz  --max-missing 0.5 --maf 0.05 --min-meanDP 10 --recode --recode-INFO-all  --stdout  >  filtered_lenient.5.vcf


NUM_VARIANTS1=$(cat filtered_lenient.5.vcf | grep -v '^#' | wc -l)
echo "Number of variants after missing and maf filters: ${NUM_VARIANTS1}"

vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" filtered_lenient.5.vcf |bgzip > filtered_lenient.6.vcf.gz

NUM_VARIANTS2=$(zcat filtered_lenient.6.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after AB  filters: ${NUM_VARIANTS2}"

#depths
vcftools --gzvcf filtered_lenient.6.vcf.gz --site-mean-depth --out filtered_lenient.6.depth

vcftools --gzvcf filtered_lenient.6.vcf.gz --max-meanDP 136 --exclude-positions ~/Tursiops-RAD-popgen/scripts/HD_exclude.txt --recode --recode-INFO-all --stdout |  bgzip > filtered_lenient.final_MitoIncluded.vcf.gz
zcat filtered_lenient.final_MitoIncluded.vcf.gz | grep -v "^#" | wc -l

# remove the mitochondrial chromosome

vcftools --gzvcf filtered_lenient.final_MitoIncluded.vcf.gz --not-chr NC_012059.1 --not-chr NW_022983433.1 --recode --recode-INFO-all --stdout |  bgzip > filtered_lenient.final.vcf.gz

bcftools annotate --set-id '%CHROM\_%POS' filtered_lenient.final.vcf.gz  -Oz -o filtered_lenient.final_ids.vcf.gz

bcftools index -t filtered_lenient.final_MitoIncluded.vcf.gz
bcftools index -t filtered_lenient.final.vcf.gz
bcftools index -t filtered_lenient.final_ids.vcf.gz


plink --vcf filtered_lenient.final_ids.vcf.gz \
	        --indep-pairwise 50 5 0.2 --allow-extra-chr --double-id \
		        --out variants_pruned

# make thinned vcf
vcftools --gzvcf filtered_lenient.final_ids.vcf.gz --snps variants_pruned.prune.in  --recode --recode-INFO-all --stdout | \
	        bgzip > filtered_lenient.final_ids_LDthin.vcf.gz

bcftools index -t filtered_lenient.final_ids_LDthin.vcf.gz

zcat filtered_lenient.final_ids_LDthin.vcf.gz | grep -v "^#" | wc -l
# 7557


