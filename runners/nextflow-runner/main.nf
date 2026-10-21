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
    if ((params.mode == 'braker3_protein' || params.mode == 'braker3') && params.orthodb && !file(params.orthodb).exists()) {
        error "--orthodb file not found: ${params.orthodb}"
    }
    if (params.mode == 'galba' && !params.prot_seq) {
        error "GALBA mode requires --prot_seq (close-relative proteomes)"
    }

    // ── Input channels ──────────────────────────────────────────────────────

    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)

    // Optional RNA-seq reads
    if (params.rnaseq) {
        reads_ch = Channel.fromFilePairs(params.rnaseq, checkIfExists: true)
    }

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
                    params.rnaseq ? reads_ch : Channel.empty(),
                    Channel.fromPath(params.orthodb, checkIfExists: true),
                    params.species,
                    'full')
            gene_gff_ch    = BRAKER3.out.gff
            gene_prot_ch   = BRAKER3.out.proteins
            gene_cds_ch    = BRAKER3.out.codingseq
    } else if (params.mode == 'braker3_protein') {
            BRAKER3(masked_genome_ch,
                    Channel.empty(),
                    Channel.fromPath(params.orthodb, checkIfExists: true),
                    params.species,
                    'protein_only')
            gene_gff_ch    = BRAKER3.out.gff
            gene_prot_ch   = BRAKER3.out.proteins
            gene_cds_ch    = BRAKER3.out.codingseq
    } else if (params.mode == 'braker3_rnaseq') {
            BRAKER3(masked_genome_ch,
                    reads_ch,
                    Channel.empty(),
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
                              params.rnaseq ? reads_ch : Channel.empty())
            FUNANNOTATE_PREDICT(FUNANNOTATE_TRAIN.out.outdir,
                                params.species,
                                params.rnaseq ? reads_ch : Channel.empty(),
                                params.prot_seq ? Channel.fromPath(params.prot_seq) : Channel.empty())
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
        EGGNOG_MAPPER(gene_prot_ch, params.eggnog_db ?: 'db/')
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