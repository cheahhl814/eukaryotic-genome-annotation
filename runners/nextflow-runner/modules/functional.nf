// ─── Module: Functional Annotation (eggNOG + InterProScan + KofamScan) ────
//
// Phase 4 of the genome annotation pipeline.
// See ../functional-annotation/SKILL.md for tiering, pitfalls, and reconciliation.

process EGGNOG_MAPPER {
    tag "${proteins.baseName}"
    label 'process_long'
    container 'quay.io/biocontainers/eggnog-mapper:2.1.15--pyh7e72e48_0'

    input:
    path proteins
    path eggnog_db                           // path to eggNOG database dir or 'db/'

    output:
    path "eggnog_annotations.tsv", emit: annotations


    script:
    """
    emapper.py \\
        -i ${proteins} \\
        --itype proteins \\
        -m diamond \\
        --tax_scope auto \\
        --data_dir ${eggnog_db} \\
        --cpu ${task.cpus} \\
        -o eggnog \\
        --output_dir eggnog_out

    mv eggnog_out/eggnog.emapper.annotations eggnog_annotations.tsv
    """
    stub:
    """
touch eggnog_annotations.tsv
    """
}

process INTERPROSCAN {
    tag "${proteins.baseName}"
    label 'process_long'
    container 'quay.io/biocontainers/interproscan:5.69_0--h7d875b9_0'

    input:
    path proteins

    output:
    path "interproscan.tsv",  emit: tsv
    path "interproscan.gff3", emit: gff3


    script:
    """
    interproscan.sh \\
        -i ${proteins} \\
        -f TSV,GFF3 \\
        -goterms \\
        -pa \\
        -cpu ${task.cpus} \\
        -o interproscan.tsv \\
        -gff3 interproscan.gff3
    """
    stub:
    """
touch interproscan.tsv interproscan.gff3
    """
}

process KOFAMSCAN {
    tag "${proteins.baseName}"
    label 'process_medium'
    container 'quay.io/biocontainers/kofamscan:1.3.0--h77072ba_2'

    input:
    path proteins

    output:
    path "kofam_annotations.tsv", emit: annotations


    script:
    """
    exec_annotation \\
        -f detail-tsv \\
        -o kofam_annotations.tsv \\
        --cpu ${task.cpus} \\
        ${proteins}
    """
    stub:
    """
touch kofam_annotations.tsv
    """
}