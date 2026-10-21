// ─── Module: Funannotate Gene Prediction (Fungal Standard) ──────────────────
//
// Phase 2 of the genome annotation pipeline — four-step fungal pipeline.
// See ../funannotate-pipeline/SKILL.md for fungal-specific pitfalls.

process FUNANNOTATE_MASK {
    tag "${genome.baseName}"
    label 'process_medium'
    container 'nextgenusfs/funannotate:latest'

    input:
    path genome

    output:
    path "masked.fa", emit: masked_genome


    script:
    """
    funannotate mask -i ${genome} -o masked.fa --cpus ${task.cpus}
    """
    stub:
    """
touch masked.fa
    """
}

process FUNANNOTATE_TRAIN {
    tag 'funannotate-train'
    label 'process_high'
    container 'nextgenusfs/funannotate:latest'

    input:
    path masked_genome
    path reads                             // paired-end reads or empty

    output:
    path "funannotate_train", emit: outdir


    script:
    def reads_flag = reads ? "-l ${reads[0]} -r ${reads[1]}" : ''

    """
    funannotate train \\
        -i ${masked_genome} \\
        -o funannotate_train \\
        ${reads_flag} \\
        --species "${params.species}" \\
        --cpus ${task.cpus}
    """
    stub:
    """
mkdir -p funannotate_train
    """
}

process FUNANNOTATE_PREDICT {
    tag "${species}"
    label 'process_long'
    container 'nextgenusfs/funannotate:latest'

    input:
    path train_outdir
    val  species
    path reads
    path prot_seq

    output:
    path "funannotate_predict", emit: outdir


    script:
    def trans_flag  = file("${train_outdir}/transcripts.fa").exists() ? "--transcript_evidence ${train_outdir}/transcripts.fa" : ''
    def prot_flag   = prot_seq ? "--protein_evidence ${prot_seq}" : ''

    """
    funannotate predict \\
        -i ${train_outdir}/masked.fa \\
        -o funannotate_predict \\
        -s "${species}" \\
        ${trans_flag} \\
        ${prot_flag} \\
        --cpus ${task.cpus}
    """
    stub:
    """
mkdir -p funannotate_predict
ln -s ../funannotate_train/masked.fa funannotate_predict/masked.fa 2>/dev/null || true
    """
}

process FUNANNOTATE_UPDATE {
    tag 'funannotate-update'
    label 'process_medium'
    container 'nextgenusfs/funannotate:latest'

    input:
    path predict_outdir

    output:
    path "funannotate_update/*.gff3",  emit: gff
    path "funannotate_update/*.faa",    emit: proteins
    path "funannotate_update/*.cds",    emit: codingseq
    path "funannotate_update",          emit: outdir


    script:
    """
    funannotate update \\
        -i ${predict_outdir} \\
        --cpus ${task.cpus}
    """
    stub:
    """
mkdir -p funannotate_update
touch funannotate_update/updated.gff3 funannotate_update/updated.faa funannotate_update/updated.cds
    """
}