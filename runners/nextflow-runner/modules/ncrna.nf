// ─── Module: ncRNA Annotation (Infernal + tRNAscan-SE + barrnap) ───────────
//
// Phase 3 of the genome annotation pipeline.
// See ../ncrna-annotation/SKILL.md for domain mode selection and pitfalls.

process INFERNAL_CMSCAN {
    tag "${genome.baseName}"
    label 'process_high'
    container 'quay.io/biocontainers/infernal:1.1.5--h77072ba_1'

    input:
    path genome

    output:
    path "ncrna.cmscan",       emit: cmscan_output
    path "ncrna.tblout",       emit: tblout
    path "ncrna.deoverlapped", emit: deoverlapped


    script:
    // Rfam database: auto-cmpress if needed
    def rfam_cm = params.rfam_cm ?: 'Rfam.cm'
    def rfam_clanin = params.rfam_clanin ?: 'Rfam.clanin'
    """
    # Calculate database size for run-comparable E-values
    RESIDUES=\$(grep -v '>' ${genome} | tr -d '\\n' | wc -c)
    ZVAL=\$(echo "scale=2; \${RESIDUES} * 2 / 1000000" | bc)

    # Ensure Rfam database is cmpressed (required for --cut_ga)
    if [ ! -f ${rfam_cm}.i1f ]; then
        cmpress ${rfam_cm}
    fi

    # Full canonical command with all correctness flags
    cmscan -Z \${ZVAL} \\
        --cut_ga --rfam --nohmmonly \\
        --tblout ncrna.tblout --fmt 2 \\
        --clanin ${rfam_clanin} \\
        ${rfam_cm} \\
        ${genome} > ncrna.cmscan

    # Mandatory deoverlap — raw table double-counts within clans
    grep -v ' = ' ncrna.tblout > ncrna.deoverlapped
    """
    stub:
    """
touch ncrna.cmscan ncrna.tblout ncrna.deoverlapped
    """
}

process TRNASCAN_SE {
    tag "${genome.baseName}"
    label 'process_low'
    container 'quay.io/biocontainers/trnascan-se:2.0.12--h7d875b9_0'

    input:
    path genome
    val  domain                              // euk | bac | arc | G

    output:
    path "trna.out",    emit: trna_output
    path "trna.gff3",   emit: trna_gff
    path "trna.stats",  emit: trna_stats


    script:
    // Domain flag: -E (euk) | -B (bac) | -A (arc) | -G (general/metagenome)
    def domain_flag = ['euk': '-E', 'bac': '-B', 'arc': '-A', 'G': '-G'].get(domain.toLowerCase().take(3), '-E')

    """
    tRNAscan-SE ${domain_flag} \\
        -o trna.out \\
        -f trna.ss \\
        -m trna.stats \\
        --gff trna.gff3 \\
        ${genome}
    """
    stub:
    """
touch trna.out trna.stats trna.gff3
    """
}

process BARRNAP {
    tag "${genome.baseName}"
    label 'process_single'
    container 'quay.io/biocontainers/barrnap:0.9--hec16e2b_0'

    input:
    path genome
    val  domain                              // euk | bac | arc | mito

    output:
    path "rrna.gff3", emit: rrna_gff


    script:
    """
    barrnap --kingdom ${domain} ${genome} > rrna.gff3
    """
    stub:
    """
touch rrna.gff3
    """
}