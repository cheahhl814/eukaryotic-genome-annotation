// ─── Module: RNA-seq Alignment (HISAT2 → samtools sort) ─────────────────────
//
// Phase 2 evidence prep — feeds BRAKER3's --bam (which requires a SORTED BAM;
// raw FASTQs are not accepted). main.nf runs this when --rnaseq is given and
// --rnaseq_bam is absent. See ../braker3-pipeline/SKILL.md and
// ../gene-prediction-decision-tree/SKILL.md.

process HISAT2_ALIGN {
    tag "${id}"
    label 'process_medium'
    // Tag verified against quay.io/biocontainers (2026-10-21)
    container 'quay.io/biocontainers/hisat2:2.2.1--h503566f_8'

    input:
    tuple val(id), path(reads)      // reads = [reads_1.fq.gz, reads_2.fq.gz]

    output:
    tuple val(id), path("aligned.sam"), emit: aligned

    script:
    """
    hisat2 -p ${task.cpus} \\
        -x ${params.hisat2_index_prefix} \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -S aligned.sam
    """
}

process SAMTOOLS_SORT {
    tag "${id}"
    label 'process_medium'
    // Tag verified against quay.io/biocontainers (2026-10-21)
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'

    input:
    tuple val(id), path(sam)

    output:
    tuple val(id), path("aligned.sorted.bam"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o aligned.sorted.bam ${sam}
    samtools index aligned.sorted.bam
    """
}