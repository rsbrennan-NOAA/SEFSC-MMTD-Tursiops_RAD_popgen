#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.indir = null
params.outdir = null
params.refDir = "/home/rbrennan/reference_genomes/bottlenose_dolphin"
params.refGenome = "${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz"
params.pattern = '*_{1,2}.fastq.gz'

// Parameter validation
if (params.indir == null) {
    error "Input directory path not specified. Please provide with --indir"
}
if (params.outdir == null) {
    error "Output directory path not specified. Please provide with --outdir"
}
//if (params.refDir == null) {
//    error "Reference directory path not specified. Please provide with --refDir"
//}
if (params.refGenome == null) {
    error "Reference genome path not specified. Please provide with --refGenome"
}

def checkRefFiles(refFiles) {
    refFiles.each { file ->
        if (!file.exists()) {
            error "Required reference file not found: ${file}. Please ensure all BWA index files exist."
        }
    }
}


// Process to align reads
process bwamem2 {
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sampleId), path(read1), path(read2)
    path(refGenome)
    path(ref_0123)
    path(ref_amb)
    path(ref_ann)
    path(ref_bwt)
    path(ref_pac)

    output:
    path("${sampleId}.bam"), emit: bam
    path("${sampleId}.bam.bai"), emit: bai

    script:
    """
    # Load modules
    module load aligners/bwa-mem2/2.2.1
    module load bio/samtools/1.19

    # Set read group
    RG="@RG\\tID:${sampleId}\\tSM:${sampleId}"

    # Align paired-end reads
    bwa-mem2 mem -t ${task.cpus} -R "\$RG" $refGenome $read1 $read2 | \
    samtools view -S -h -u - | \
    samtools sort -  > ${sampleId}.bam

    # Index BAM file
    samtools index ${sampleId}.bam
    """
}

// Workflow
workflow {
    // Input channel with sample ID extraction for paired-end reads
    reads_ch = channel
        .fromFilePairs("${params.indir}/${params.pattern}")
        .map { sampleId, files -> 
            return tuple(sampleId, files[0], files[1])
        }

    // Reference genome and associated files
    refGenome = file(params.refGenome)
    ref_0123 = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.0123")
    ref_amb = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.amb")
    ref_ann = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.ann")
    ref_bwt = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.bwt.2bit.64")
    ref_pac = file("${params.refDir}/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz.pac")
    
    // Check if all reference files exist
    checkRefFiles([refGenome, ref_0123, ref_amb, ref_ann, ref_bwt, ref_pac])

    // Run alignment
    bwamem2(reads_ch, refGenome, ref_0123, ref_amb, ref_ann, ref_bwt, ref_pac)
}
