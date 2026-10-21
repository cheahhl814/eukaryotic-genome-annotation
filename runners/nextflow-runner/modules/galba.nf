// ─── Module: GALBA Gene Prediction ───────────────────────────────────────────
//
// Phase 2 of the genome annotation pipeline — protein-only with close relatives.
// See ../galba-pipeline/SKILL.md for when to prefer GALBA over BRAKER3.

process GALBA {
    tag "${species}"
    label 'process_long'
    container 'quay.io/biocontainers/galba:1.0.3--h9ee8adc_0'

    input:
    path genome
    path prot_seq
    val  species

    output:
    path "galba_out/galba.gff3",       emit: gff
    path "galba_out/galba.aa",          emit: proteins
    path "galba_out/galba.codingseq",   emit: codingseq
    path "galba_out",                   emit: outdir


    script:
    """
    galba.pl \\
        --genome=${genome} \\
        --prot_seq=${prot_seq} \\
        --species=${species} \\
        --threads=${task.cpus} \\
        --workingdir=galba_out
    """
    stub:
    """
mkdir -p galba_out
touch galba_out/galba.gff3 galba_out/galba.aa galba_out/galba.codingseq
    """
}