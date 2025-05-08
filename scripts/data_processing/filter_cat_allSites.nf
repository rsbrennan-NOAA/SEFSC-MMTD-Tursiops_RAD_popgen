
// Define parameters
params.outdir = "${HOME}/Tursiops-RAD-popgen/analysis/variants/"
params.indir = "${HOME}/Tursiops-RAD-popgen/analysis/variants/all_sites"

process CONCAT_VCF {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path vcf_dir    
    
    output:
    path 'variants_allSites_raw_merged.vcf.gz', emit: merged_vcf
    path 'variants_allSites_raw_merged.vcf.gz.tbi'
    
    script:
    """
    module load bio/bcftools/1.11
    module load bio/htslib/1.19
    ls ${vcf_dir}/variants_raw_allSites*.vcf.gz > vcf_allSites_list.txt

    bcftools concat -f vcf_allSites_list.txt -O z -o variants_allSites_raw_merged.vcf.gz

    tabix -p vcf variants_allSites_raw_merged.vcf.gz
    """
}

process FILTER_VCF {
        publishDir params.outdir, mode: 'copy'

        input:
        path merged_vcf

        output:
        path 'filtered.invariant.vcf.gz', emit: filtered_invariant

        script:
        """
	source /opt/bioinformatics/mambaforge/etc/profile.d/conda.sh

	conda activate vcflib-1.0.9
	module load bio/vcftools
	module load bio/bcftools
        
	vcftools --gzvcf ${merged_vcf} \
        --max-maf 0.05 \
        --remove-indels \
        --remove ~/Tursiops-RAD-popgen/scripts/rm_missing.txt \
	--max-missing 0.7 \
        --min-meanDP 10 \
        --max-meanDP 136 \
	--not-chr NC_012059.1 \
        --recode --recode-INFO-all --stdout | bgzip -c > filtered.invariant.vcf.gz

        """
}

process COMBINE_VCFS {
    publishDir params.outdir, mode: 'move'
    input:
    path filtered_invariant
    output:
    path 'combined_filtered_invariant.vcf.gz', emit: combined_filtered
    path 'combined_filtered_invariant.vcf.gz.tbi'

    script:
    """
    source /opt/bioinformatics/mambaforge/etc/profile.d/conda.sh
    conda activate vcflib-1.0.9

    module load bio/bcftools
    module load bio/htslib/1.19
    module load bio/vcftools

    final_vcf="${HOME}/Tursiops-RAD-popgen/analysis/variants/filtered.final.vcf.gz"
    
    tabix -f ${filtered_invariant}
    tabix -f \${final_vcf}

   # Combine the two VCFs
    bcftools concat \\
        --allow-overlaps \\
        \${final_vcf} ${filtered_invariant} \\
        -O z -o combined_filtered_invariant.1.vcf.gz
    
    vcftools --gzvcf combined_filtered_invariant.1.vcf.gz \
	--chr NC_047034.1 \
	--chr NC_047035.1 \
	--chr NC_047036.1 \
	--chr NC_047037.1 \
	--chr NC_047038.1 \
	--chr NC_047039.1 \
	--chr NC_047040.1 \
	--chr NC_047041.1 \
	--chr NC_047042.1 \
	--chr NC_047043.1 \
	--chr NC_047044.1 \
	--chr NC_047045.1 \
	--chr NC_047046.1 \
	--chr NC_047047.1 \
	--chr NC_047048.1 \
	--chr NC_047049.1 \
	--chr NC_047050.1 \
	--chr NC_047051.1 \
	--chr NC_047052.1 \
	--chr NC_047053.1 \
	--chr NC_047054.1 \
	--chr NC_047055.1 \
	--recode --recode-INFO-all --stdout | bgzip -c > combined_filtered_invariant.vcf.gz

    tabix -p vcf combined_filtered_invariant.vcf.gz
    """
}



// Main workflow
workflow {
    vcf_dir = Channel.fromPath(params.indir, checkIfExists: true)
    CONCAT_VCF(vcf_dir)
    FILTER_VCF(CONCAT_VCF.out.merged_vcf)
    COMBINE_VCFS(FILTER_VCF.out.filtered_invariant)
}


