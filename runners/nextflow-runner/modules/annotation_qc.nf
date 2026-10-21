// ─── Module: Annotation QC (BUSCO + compleasm + OMArk + sanity panel) ─────
//
// Phase 5 of the genome annotation pipeline.
// See ../annotation-qc/SKILL.md for interpretation logic and pitfalls.
// The workflow produces the metrics; the agent interprets them.

process BUSCO_PROTEOME {
    tag "${lineage}"
    label 'process_medium'
    container 'quay.io/biocontainers/busco:5.8.0--h7d875b9_0'

    input:
    path proteins
    val  lineage

    output:
    path "busco_prot", emit: outdir


    script:
    """
    busco -i ${proteins} -m proteins -l ${lineage} -o busco_prot -c ${task.cpus}
    """
    stub:
    """
mkdir -p busco_prot
    """
}

process BUSCO_GENOME {
    tag "${lineage}"
    label 'process_huge'
    container 'quay.io/biocontainers/busco:5.8.0--h7d875b9_0'

    input:
    path genome
    val  lineage

    output:
    path "busco_genome", emit: outdir


    script:
    """
    busco -i ${genome} -m genome -l ${lineage} -o busco_genome -c ${task.cpus}
    """
    stub:
    """
mkdir -p busco_genome
    """
}

process COMPLEASM {
    tag "${lineage}"
    label 'process_medium'
    container 'quay.io/biocontainers/compleasm:0.2.7--h77072ba_0'

    input:
    path genome
    path proteins
    val  lineage

    output:
    path "compleasm_out", emit: outdir


    script:
    """
    # Genome mode (assembly capacity)
    compleasm run -a ${genome} -l ${lineage} -o compleasm_genome -t ${task.cpus}

    # Protein mode (delivered proteome)
    compleasm protein -p ${proteins} -l ${lineage} -o compleasm_protein -t ${task.cpus}
    """
    stub:
    """
mkdir -p compleasm_out
    """
}

process OMARK {
    tag "omark"
    label 'process_long'
    container 'quay.io/biocontainers/omark:1.1.0--h7d875b9_0'

    input:
    path proteins

    output:
    path "omark_out", emit: outdir


    script:
    """
    omamer search --db \${OMARK_DB:-LUCA.h5} --query ${proteins} --out proteins.omamer
    omark -f proteins.omamer -d \${OMARK_DB:-LUCA.h5} -o omark_out
    """
    stub:
    """
mkdir -p omark_out
    """
}

process SANITY_PANEL {
    tag "sanity"
    label 'process_single'
    container 'quay.io/biocontainers/gffutils:0.12--pyh7e72e48_0'

    input:
    path gff

    output:
    path "sanity_panel.tsv", emit: metrics


    script:
    """
    python3 << 'EOF'
    import gffutils
    import csv

    db = gffutils.create_db("${gff}", ':memory:', merge_strategy='merge')
    genes = list(db.features_of_type('gene'))
    mrnas = list(db.features_of_type(['mRNA', 'transcript']))
    exon_counts = [len(list(db.children(tx, featuretype='exon'))) for tx in mrnas]

    mono_frac = sum(1 for e in exon_counts if e == 1) / len(exon_counts) if exon_counts else 0
    mrna_per_gene = len(mrnas) / len(genes) if genes else 0

    with open("sanity_panel.tsv", "w") as f:
        w = csv.writer(f, delimiter="\\t")
        w.writerow(["metric", "value"])
        w.writerow(["gene_count", len(genes)])
        w.writerow(["mrna_per_gene", f"{mrna_per_gene:.4f}"])
        w.writerow(["mono_exonic_fraction", f"{mono_frac:.4f}"])
        if mrna_per_gene <= 1.001:
            w.writerow(["warning", "mRNA:gene == 1.00 -- isoform/UTR-naive"])
        if mono_frac > 0.30:
            w.writerow(["warning", f"mono-exonic fraction {mono_frac:.1%} high -- check masking/contamination"])
    EOF
    """
    stub:
    """
touch sanity_panel.tsv
    """
}