#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Input params
params.vcf_file = "${HOME}/Tursiops-RAD-popgen/analysis/variants/filtered.final.noMAF.vcf.gz"
params.num_bootstraps = 100
params.output_dir = "${HOME}/Tursiops-RAD-popgen/analysis/moments/bootstrap_sfs"
params.genome_file = "${HOME}/reference_genomes/bottlenose_dolphin/GCF_011762595.1_mTurTru1.mat.Y_genomic_genomeFile.txt"
params.pop_file = "${HOME}/Tursiops-RAD-popgen/analysis/pop_structure/fourpop_all_noHybrids.clust"
params.proj = "66,64,66,66"
params.pop_order = "Coastal_Atlantic,Coastal_Gulf,Intermediate,Offshore"
params.chunks_dir = "${HOME}/Tursiops-RAD-popgen/analysis/moments/chunks"
params.num_parts = 5

// Print parameters to console
log.info """\
        BOOTSTRAP VCF for moments
        --------------------------------------------
        vcf_file        : ${params.vcf_file}
        num_bootstraps  : ${params.num_bootstraps}
        output_dir      : ${params.output_dir}
        genome_file     : ${params.genome_file}
        chunks_dir      : ${params.chunks_dir}
        num_parts       : ${params.num_parts}
        pop_file        : ${params.pop_file}
        proj            : ${params.proj}
        pop_order       : ${params.pop_order}
        """
        .stripIndent()

process PROCESSVCF {
    tag "process_vcf"
    publishDir params.chunks_dir, mode: 'copy', pattern: '*_part*.{bed,vcf}'
    
    input:
    path vcf
    path genome
    
    output:
    path 'header.txt', emit: header
    path '*_part*.vcf', emit: chunks
    
    script:
    """

    module load bio/bedtools/2.31.1
    module load bio/vcftools/0.1.16

    # Extract header from VCF
    zcat ${vcf} | grep "^#" > header.txt
    
    # Extract chromosomes from VCF
    zcat ${vcf} | grep -v "^#" | cut -f1 | sort | uniq > chromosomes.txt
    
    # Process each chromosome
    while read chromosome; do
        # Get chromosome length
        CHR_LENGTH=\$(grep -w "\$chromosome" ${genome} | cut -f2)
        
        # Calculate chunk sizes
        PART_SIZE=\$((CHR_LENGTH / ${params.num_parts}))
        
        # Create modified chromosome name
        modChr=\$(echo \$chromosome | tr '.' '_')
        
  # Create chunks
        for PART in \$(seq 0 \$((${params.num_parts}-1))); do
            START=\$((PART * PART_SIZE))
            # For the last part, ensure we go to the end of the chromosome
            if [ \$PART -eq \$((${params.num_parts}-1)) ]; then
                END=\$CHR_LENGTH
            else
                END=\$(((PART + 1) * PART_SIZE))
            fi
            
            # make chunk ID
            CHUNK_ID="\${modChr}_part\$((PART + 1))"
            
            # make bed file for this chunk
            echo -e "\$chromosome\\t\${START}\\t\${END}" > \${CHUNK_ID}.bed
            
            # Extract VCF for this chunk
            bedtools intersect -a ${vcf} -b \${CHUNK_ID}.bed -header > \${CHUNK_ID}.vcf
        done
    done < chromosomes.txt
    """
}

// Generate bootstrap replicates
process GENBOOTSTRAP {
    tag "bootstrap_${bootstrap_id}" //  each bootstrap replicate gets labeled with its number 
    publishDir "${params.output_dir}/vcfs", mode: 'copy'
    
    input:
    path chunks
    path header
    each bootstrap_id
    
    output:
    tuple val(bootstrap_id), path("bootstrap_${bootstrap_id}.vcf"), emit: bootstrap_vcf
    
    script:
    """
    # Start with header, then will append variants after
    cat ${header} > bootstrap_${bootstrap_id}.vcf
    
    # Store chunk files in an array using the chunks variable
    chunk_files=(${chunks})
    num_chunks=\${#chunk_files[@]}
    
    # Sample chunks with replacement
    for j in \$(seq 1 \$num_chunks); do
        # Select a random chunk
        random_index=\$((RANDOM % num_chunks))
        random_chunk=\${chunk_files[\$random_index]}
        
        # Append to bootstrap file
        grep -v "^#" \$random_chunk >> bootstrap_${bootstrap_id}.vcf
    done
    """
}

// Run easySFS on each bootstrap replicate
process EASYSFS {
    tag "easySFS_${bootstrap_id}"
    publishDir "${params.output_dir}/bootstrap_${bootstrap_id}", mode: 'copy'
    
    input:
    tuple val(bootstrap_id), path(bootstrap_vcf)
    path pop_file
    
    output:
    path "*"
    
    script:
    """
    source /opt/bioinformatics/mambaforge/etc/profile.d/conda.sh
    conda activate easySFS
    
    # Run easySFS
    ~/bin/easySFS/easySFS.py -i ${bootstrap_vcf} -p ${pop_file} \
        --order ${params.pop_order} --proj ${params.proj} -a -f -o .
    """
}
workflow {
    // Create input channels
    vcf_ch = Channel.fromPath(params.vcf_file)
                    .ifEmpty { error "Cannot find VCF file: ${params.vcf_file}" }
    
    genome_ch = Channel.fromPath(params.genome_file)
                      .ifEmpty { error "Cannot find genome file: ${params.genome_file}" }
    
    pop_ch = Channel.fromPath(params.pop_file)
                   .ifEmpty { error "Cannot find population file: ${params.pop_file}" }
    
    // Create bootstrap IDs channel
    bootstrap_ids_ch = Channel.of(1..params.num_bootstraps)
    
    // Process VCF to get header and chunks
    processvcf_results = PROCESSVCF(vcf_ch, genome_ch)
   
    // Process channels for bootstrap generation
    chunks_processed = processvcf_results.chunks
                            .flatten()
                            .filter { it.toString().endsWith('.vcf') }
                            .collect()

    // Generate bootstrap replicates
   bootstrap_results = GENBOOTSTRAP(
    chunks_processed, 
    processvcf_results.header, 
    bootstrap_ids_ch
) 
    // Run easySFS on each bootstrap replicate
    EASYSFS(bootstrap_results.bootstrap_vcf, pop_ch.collect())
}
