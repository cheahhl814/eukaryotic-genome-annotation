---
name: "eukaryotic-genome-annotation/nextflow-runner"
description: "Nextflow DSL2 workflow for the genome annotation pipeline. Handles reproducible, containerized execution of repeat masking, gene prediction, ncRNA annotation, functional annotation, and QC. The agentic skills handle routing, interpretation, and failure recovery; this workflow handles unattended execution."
version: 2
triggers:
  - "run annotation workflow"
  - "nextflow annotate genome"
created: "2026-07-11"
updated: "2026-10-21"
---

# Genome Annotation Pipeline — Nextflow Workflow

**"Run the genome annotation pipeline unattended"** → Use this workflow for reproducible, containerized execution of the annotation pipeline. The agentic skills (parent `SKILL.md` and sub-skills) handle routing, QC interpretation, and failure recovery; this workflow handles deterministic execution with checkpointing and resume.

## When to Use the Workflow vs Manual Execution

| Scenario                                          | Use workflow | Use manual/agent |
| ------------------------------------------------- | ------------ | ---------------- |
| Production annotation run on a cluster            | ✅            |                  |
| Need checkpointing/resume for long runs           | ✅            |                  |
| Re-annotating multiple genomes with same settings | ✅            |                  |
| First-time annotation, need to interpret QC gates |              | ✅                |
| Assembly might not pass Phase 0 gate              |              | ✅                |
| Deciding which predictor to use                   |              | ✅                |
| Diagnosing a BRAKER failure mid-run               |              | ✅                |
| Curating a repeat library (requires judgment)     |              | ✅                |
| Interpreting "is 40% Unknown acceptable?"         |              | ✅                |

**Pattern:** Use the agentic skill for Phase 0 (assembly QC gate), Phase 1 (repeat library curation), and Phase 5 (QC interpretation). Use the workflow for the mechanical execution of Phases 1–4 and Phase 5's BUSCO/OMArk runs.

## Quick Start

```bash
# Phase 0: Assembly QC gate (run manually/agent FIRST)
busco -i assembly.fa -m genome -l vertebrata_odb12 -o busco_asm -c 16

# If assembly passes, run the workflow:
nextflow run main.nf \
  -profile singularity,slurm \
  --genome assembly.fa \
  --mode braker3 \
  --species 'My_species' \
  --orthodb OrthoDB/Vertebrata.fa \
  --rnaseq 'reads/*_{1,2}.fq.gz' \
  --outdir results/

# Resume from checkpoint if interrupted:
nextflow run main.nf -resume -profile singularity,slurm
```

## Pipeline Modes

| `--mode`          | Gene predictor                   | Required inputs                     |
| ----------------- | -------------------------------- | ----------------------------------- |
| `braker3`         | BRAKER3 full (RNA-seq + protein) | `--genome`, `--rnaseq`+`--hisat2_index_prefix` (or `--rnaseq_bam`), `--orthodb` (optional `--skip_ncrna`/`--skip_functional` with their DBs) |
| `braker3_protein` | BRAKER3 protein-only             | `--genome`, `--orthodb`             |
| `braker3_rnaseq`  | BRAKER3 RNA-seq-only             | `--genome`, `--rnaseq`+`--hisat2_index_prefix`, or `--rnaseq_bam` |
| `galba`           | GALBA (close-relative proteins)  | `--genome`, `--prot_seq`            |
| `funannotate`     | Funannotate (fungal)             | `--genome`, `--rnaseq` (optional)   |

> [!NOTE]
> **BRAKER3 takes a sorted BAM, not FASTQs.** Raw RNA-seq FASTQs are aligned through the
> built-in `HISAT2_ALIGN → SAMTOOLS_SORT` stage (`modules/rnaseq_align.nf`; containers
> verified on quay.io biocontainers) before BRAKER3's `--bam`. Have a pre-aligned BAM?
> Pass `--rnaseq_bam` and the alignment stage is skipped. The GALBA protein-evidence route
> can borrow pre-made reference proteomes from the biodb-fetch skill (see the skill root's
> `galba-pipeline/SKILL.md` §Prerequisites).

## Phase Control

```bash
# Skip phases if already completed:
--skip_repeat        # Use pre-masked genome (--genome_masked)
--skip_ncrna         # Skip ncRNA annotation
--skip_functional    # Skip eggNOG/InterProScan/KofamScan
--skip_qc            # Skip final QC (run manually)

# Provide pre-computed inputs:
--genome_masked masked.fa     # Skip repeat masking
--braker_gff braker.gff3     # Skip gene prediction
--proteins proteins.faa       # Skip gene prediction for functional annotation
```

## Output Structure

```
results/
├── repeat_annotation/
│   ├── repeatmasker/
│   │   ├── genome.fa.masked          # Soft-masked genome
│   │   ├── genome.fa.out             # Repeat table
│   │   └── genome.fa.tbl             # Summary statistics
│   └── repeatmodeler/
│       └── mydb-families.fa          # De novo repeat library
├── gene_prediction/
│   ├── braker/                       # or galba/ or funannotate/
│   ├── braker.gff3                   # Gene models
│   ├── braker.aa                     # Predicted proteins
│   └── braker.codingseq              # Coding sequences
├── ncrna/
│   ├── infernal/
│   │   └── ncrna.cmscan.tblout       # Structured ncRNA hits
│   ├── trna/
│   │   ├── trna.gff3                 # tRNA annotations
│   │   └── trna.stats                # tRNA statistics
│   └── rrna/
│       └── rrna.gff3                 # rRNA annotations
├── functional/
│   ├── eggnog/
│   │   └── eggnog_annotations.tsv    # eggNOG-mapper output
│   ├── interproscan/
│   │   └── interproscan.tsv          # InterProScan output
│   └── kofam/
│       └── kofam_annotations.tsv     # KofamScan output
├── qc/
│   ├── busco_proteome/               # BUSCO on predicted proteins
│   ├── busco_genome/                 # BUSCO on assembly
│   ├── compleasm/                    # compleasm output
│   ├── omark/                        # OMArk output
│   └── sanity_panel.tsv              # Gene-set sanity metrics
└── pipeline_info/                    # Nextflow trace/report/timeline
```

## Pitfalls

- **Workflow does NOT skip Phase 0:** Always run assembly QC manually before the workflow. The workflow includes a BUSCO genome check but does not gate on it — the agent must decide.
- **Repeat library curation requires judgment:** The workflow runs RepeatModeler2 + RepeatMasker, but library decontamination (removing domesticated genes) requires manual BLAST review. Use `--skip_repeat` with a curated library if you need this step.
- **BRAKER3 failures need manual intervention:** GeneMark key expiration, special characters in FASTA headers, and wrong OrthoDB partition all require agent judgment. The workflow will fail and log the error.
- **InterProScan is slow:** On large proteomes (>50k proteins), InterProScan can take days. Consider chunking or running with `-dp` for local compute.

## Verification

- [ ] Assembly passed Phase 0 QC before running the workflow.
- [ ] Correct `--mode` selected by decision tree.
- [ ] OrthoDB partition is the smallest clade containing target species.
- [ ] `--genome_masked` provided if skipping repeat masking.
- [ ] Workflow completed without errors (`nextflow log`).
- [ ] QC outputs inspected manually (agent interprets BUSCO/OMArk/sanity panel).