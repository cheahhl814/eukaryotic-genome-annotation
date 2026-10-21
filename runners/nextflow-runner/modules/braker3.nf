// ─── Module: BRAKER3 Gene Prediction ────────────────────────────────────────
//
// Phase 2 of the genome annotation pipeline.
// Supports three modes: full (RNA-seq + protein), protein-only, RNA-seq-only.
// See ../braker3-pipeline/SKILL.md for mode selection, pitfalls, and verification.

process BRAKER3 {
    tag "${species}"
    label 'process_long'
    container 'teambraker/braker3:latest'

    input:
    path genome
    path rnaseq_bam                       // empty if no RNA-seq
    path orthodb                           // empty if no protein evidence
    val  species
    val  mode                              // 'full' | 'protein_only' | 'rnaseq_only'

    output:
    path "braker3_out/braker.gff3",      emit: gff
    path "braker3_out/braker.aa",         emit: proteins
    path "braker3_out/braker.codingseq",  emit: codingseq
    path "braker3_out",                   emit: outdir


    script:
    // Mode-specific flags
    def prot_flag  = mode != 'rnaseq_only' ? "--prot_seq ${orthodb}" : ''
    def bam_flag   = mode != 'protein_only' && rnaseq_bam ? "--bam ${rnaseq_bam}" : ''

    """
    # Strip special characters from FASTA headers (common BRAKER failure)
    sed 's/ .*//' ${genome} > genome_clean.fa

    braker.pl \\
        --genome=genome_clean.fa \\
        ${prot_flag} \\
        ${bam_flag} \\
        --softmasking \\
        --threads=${task.cpus} \\
        --species=${species} \\
        --gff3 \\
        --workingdir=braker3_out
    """
    stub:
    """
mkdir -p braker3_out
touch braker3_out/braker.gff3 braker3_out/braker.aa braker3_out/braker.codingseq
    """
}