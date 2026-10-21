#!/usr/bin/env nextflow
/*
 * Genome Annotation Pipeline
 * ===========================
 * Reproducible, containerized eukaryotic genome annotation.
 *
 * Phases:
 *   1. Repeat annotation & soft-masking
 *   2. Gene prediction (BRAKER3 / GALBA / Funannotate)
 *   3. ncRNA annotation (Infernal + tRNAscan-SE + barrnap)
 *   4. Functional annotation (eggNOG-mapper + InterProScan + KofamScan)
 *   5. Annotation QC (BUSCO + compleasm + OMArk + sanity panel)
 *
 * Agent decisions (routing, QC interpretation, library curation) are handled
 * by the agentic skills. This workflow handles deterministic execution.
 *
 * See workflow/SKILL.md for when to use workflow vs manual execution.
 */

nextflow.enable.dsl = 2

// ─── Modules ────────────────────────────────────────────────────────────────
// ─── Modules ────────────────────────────────────────────────────────────────

include { REPEATMODELER2; REPEATMASKER } from './modules/repeatmasker'
include { BRAKER3 }                          from './modules/braker3'
include { HISAT2_ALIGN; SAMTOOLS_SORT }      from './modules/rnaseq_align'
include { GALBA }                            from './modules/galba'
include { FUNANNOTATE_MASK; FUNANNOTATE_TRAIN; FUNANNOTATE_PREDICT; FUNANNOTATE_UPDATE } from './modules/funannotate'
include { INFERNAL_CMSCAN; TRNASCAN_SE; BARRNAP } from './modules/ncrna'
include { EGGNOG_MAPPER; INTERPROSCAN; KOFAMSCAN } from './modules/functional'
include { BUSCO_PROTEOME; BUSCO_GENOME; COMPLEASM; OMARK; SANITY_PANEL } from './modules/annotation_qc'

// ─── Workflow ───────────────────────────────────────────────────────────────

workflow {

    // ─── Parameter validation (must live inside the workflow — DSL2 forbids
    //     mixing top-level statements with script declarations; `error` aborts) ──

    if (!params.genome) {
        error "Missing required parameter: --genome <path/to/assembly.fa>"
    }
    if (!['braker3', 'braker3_protein', 'braker3_rnaseq', 'galba', 'funannotate'].contains(params.mode)) {
        error "Invalid --mode: ${params.mode}. Use: braker3, braker3_protein, braker3_rnaseq, galba, funannotate"
    }
    if (params.mode == 'braker3' && !params.rnaseq && !params.orthodb) {
        error "BRAKER3 full mode requires --rnaseq and/or --orthodb"
    }
    if (params.mode == 'galba' && !params.prot_seq) {
        error "GALBA mode requires --prot_seq (close-relative proteomes)"
    }
    if (params.mode == 'braker3_rnaseq' && !params.rnaseq && !params.rnaseq_bam) {
        error "BRAKER3 rnaseq_only mode requires --rnaseq (FASTQs + --hisat2_index_prefix) or --rnaseq_bam (pre-aligned)"
    }
    if (params.rnaseq && !params.rnaseq_bam && !params.hisat2_index_prefix) {
        error "FASTQ RNA-seq input requires --hisat2_index_prefix (hisat2-build output prefix). Alternatively supply --rnaseq_bam."
    }
    if (!params.skip_ncrna && !params.rfam_cm) {
        error "Phase 3 requires --rfam_cm (path to Rfam.cm; fetch Rfam.clanin too for deoverlap) or --skip_ncrna"
    }
    if (!params.skip_functional && !params.eggnog_db) {
        error "Phase 4 requires --eggnog_db (eggNOG data dir via download_eggnog_data.py) or --skip_functional"
    }
    if (params.orthodb && !file(params.orthodb).exists()) {
        error "--orthodb file not found: ${params.orthodb}"
    }
    if (params.mode == 'braker3_rnaseq' && !params.rnaseq && !params.rnaseq_bam) {
        error "BRAKER3 rnaseq_only mode requires --rnaseq (FASTQs + --hisat2_index_prefix) or --rnaseq_bam (pre-aligned)"
    }
    if (params.rnaseq && !params.rnaseq_bam && !params.hisat2_index_prefix) {
        error "FASTQ RNA-seq input requires --hisat2_index_prefix (hisat2-build output prefix). Alternatively supply --rnaseq_bam."
    }

    // ── Input channels ──────────────────────────────────────────────────────

    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)

    // ── Optional RNA-seq reads → BRAKER3 needs a SORTED BAM; funannotate takes FASTQs ──

    // Evidence packs: tuple(has_evidence, file). Never pass Channel.empty() as an
    // "optional" input — an empty input channel zeroes the cross-product and the
    // process never fires (the DAG deadlocks). Placeholders are created once, here.
    PH_DIR = file("${launchDir}/.nf-placeholders", checkIfExists: false)
    if (!PH_DIR.exists()) PH_DIR.mkdirs()
    // nf2 strict syntax: Groovy `for` is not supported in a workflow body — use .each
    ['empty.bam', 'empty.faa', 'empty.fastq'].each { ph ->
        def phf = file("${launchDir}/.nf-placeholders/${ph}", checkIfExists: false)
        if (!phf.exists()) phf.text = ''
    }
    EMPTY_BAM   = file("${launchDir}/.nf-placeholders/empty.bam")
    EMPTY_FAA   = file("${launchDir}/.nf-placeholders/empty.faa")
    EMPTY_READS = file("${launchDir}/.nf-placeholders/empty.fastq")

    if (params.rnaseq)  reads_ch = Channel.fromFilePairs(params.rnaseq, checkIfExists: true)
    reads_pack = params.rnaseq
        ? reads_ch.map { id, pair -> tuple(true, pair) }
        : Channel.value(EMPTY_READS).map { f -> tuple(false, [f]) }

    bam_pack = params.rnaseq
        ? SAMTOOLS_SORT(HISAT2_ALIGN(reads_ch).aligned).bam.map { id, pf -> tuple(true, pf) }
        : params.rnaseq_bam
            ? Channel.value(file(params.rnaseq_bam, checkIfExists: true)).map { pf -> tuple(true, pf) }
            : Channel.value(EMPTY_BAM).map { pf -> tuple(false, pf) }

    prot_pack = params.orthodb
        ? Channel.value(file(params.orthodb, checkIfExists: true)).map { pf -> tuple(true, pf) }
        : Channel.value(EMPTY_FAA).map { pf -> tuple(false, pf) }

    // ── Phase 1: Repeat annotation & soft-masking ──────────────────────────

    if (!params.skip_repeat && !params.genome_masked) {
        REPEATMODELER2(genome_ch)
        REPEATMASKER(genome_ch, REPEATMODELER2.out.library)
        masked_genome_ch = REPEATMASKER.out.masked_genome
    } else if (params.genome_masked) {
        masked_genome_ch = Channel.fromPath(params.genome_masked, checkIfExists: true)
    } else {
        masked_genome_ch = genome_ch  // use unmasked if explicitly skipped
    }

    // ── Phase 2: Gene prediction ───────────────────────────────────────────

    if (params.mode == 'braker3') {
            BRAKER3(masked_genome_ch,
                    bam_pack,
                    prot_pack,
                    params.species,
                    'full')
            gene_gff_ch    = BRAKER3.out.gff
            gene_prot_ch   = BRAKER3.out.proteins
            gene_cds_ch    = BRAKER3.out.codingseq
    } else if (params.mode == 'braker3_protein') {
            BRAKER3(masked_genome_ch,
                    Channel.value(file(EMPTY_BAM)).map { pf -> tuple(false, pf) },
                    prot_pack,
                    params.species,
                    'protein_only')
            gene_gff_ch    = BRAKER3.out.gff
            gene_prot_ch   = BRAKER3.out.proteins
            gene_cds_ch    = BRAKER3.out.codingseq
    } else if (params.mode == 'braker3_rnaseq') {
            BRAKER3(masked_genome_ch,
                    bam_pack,
                    Channel.value(file(EMPTY_FAA)).map { pf -> tuple(false, pf) },
                    params.species,
                    'rnaseq_only')
            gene_gff_ch    = BRAKER3.out.gff
            gene_prot_ch   = BRAKER3.out.proteins
            gene_cds_ch    = BRAKER3.out.codingseq
    } else if (params.mode == 'galba') {
            GALBA(masked_genome_ch,
                  Channel.fromPath(params.prot_seq, checkIfExists: true),
                  params.species)
            gene_gff_ch    = GALBA.out.gff
            gene_prot_ch   = GALBA.out.proteins
            gene_cds_ch    = GALBA.out.codingseq
    } else if (params.mode == 'funannotate') {
            FUNANNOTATE_MASK(masked_genome_ch)
            FUNANNOTATE_TRAIN(FUNANNOTATE_MASK.out.masked_genome,
                              reads_pack)
            FUNANNOTATE_PREDICT(FUNANNOTATE_TRAIN.out.outdir,
                                params.species,
                                reads_pack,
                                prot_pack)
            FUNANNOTATE_UPDATE(FUNANNOTATE_PREDICT.out.outdir)
            gene_gff_ch    = FUNANNOTATE_UPDATE.out.gff
            gene_prot_ch   = FUNANNOTATE_UPDATE.out.proteins
            gene_cds_ch    = FUNANNOTATE_UPDATE.out.codingseq
    } else {
            error "Unknown mode: ${params.mode}. Use: braker3, braker3_protein, braker3_rnaseq, galba, funannotate"
    }

    // Override with pre-computed gene models if provided
    if (params.braker_gff) {
        gene_gff_ch  = Channel.fromPath(params.braker_gff, checkIfExists: true)
        gene_prot_ch = Channel.fromPath(params.proteins, checkIfExists: true)
    }

    // ── Phase 3: ncRNA annotation ──────────────────────────────────────────

    if (!params.skip_ncrna) {
        INFERNAL_CMSCAN(genome_ch)
        TRNASCAN_SE(genome_ch, params.ncrna_domain ?: 'euk')
        BARRNAP(genome_ch, params.ncrna_domain ?: 'euk')
    }

    // ── Phase 4: Functional annotation ─────────────────────────────────────

    if (!params.skip_functional) {
        EGGNOG_MAPPER(gene_prot_ch, params.eggnog_db)
        INTERPROSCAN(gene_prot_ch)
        KOFAMSCAN(gene_prot_ch)
    }

    // ── Phase 5: Annotation QC ─────────────────────────────────────────────

    if (!params.skip_qc) {
        BUSCO_PROTEOME(gene_prot_ch, params.busco_lineage ?: 'vertebrata_odb12')
        BUSCO_GENOME(genome_ch, params.busco_lineage ?: 'vertebrata_odb12')
        COMPLEASM(genome_ch, gene_prot_ch, params.busco_lineage ?: 'vertebrata')
        OMARK(gene_prot_ch)
        SANITY_PANEL(gene_gff_ch)
    }
}