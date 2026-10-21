// ─── Module: Repeat Annotation (RepeatModeler2 + RepeatMasker) ──────────────
//
// Phase 1 of the genome annotation pipeline.
// See ../repeat-annotation/SKILL.md for decision logic and pitfalls.

process REPEATMODELER2 {
    tag "${genome.baseName}"
    label 'process_high'
    container 'quay.io/biocontainers/repeatmodeler:2.0.5--h2f7e887_0'

    input:
    path genome

    output:
    path 'repeatmodeler/*/consensi.fa.classified', emit: library
    path 'repeatmodeler', emit: outdir


    script:
    """
    BuildDatabase -name ${genome.baseName} ${genome}
    RepeatModeler -database ${genome.baseName} -threads ${task.cpus} -LTRStruct -pa ${task.cpus}
    """
    stub:
    """
mkdir -p repeatmodeler/sub
touch repeatmodeler/sub/consensi.fa.classified
    """
}

process REPEATMASKER {
    tag "${genome.baseName}"
    label 'process_high'
    container 'quay.io/biocontainers/repeatmasker:4.1.7--pl5261h570b1f4_0'

    input:
    path genome
    path library

    output:
    path "${genome.baseName}.masked", emit: masked_genome
    path "${genome.baseName}.out",     emit: repeat_table
    path "${genome.baseName}.tbl",     emit: summary_stats


    script:
    """
    # Soft-mask: lowercase repeats (MANDATORY for gene prediction)
    # -nolow: skip low-complexity/simple repeats that overlap coding microsatellites
    RepeatMasker -lib ${library} -xsmall -gff -e rmblast -pa ${task.cpus} -nolow ${genome}
    """
    stub:
    """
touch ${genome.baseName}.masked ${genome.baseName}.out ${genome.baseName}.tbl
    """
}