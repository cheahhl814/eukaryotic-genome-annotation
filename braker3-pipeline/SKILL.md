---
version: 2
updated: "2026-06-29"
triggers:
  - "braker3"
  - "BRAKER3 pipeline"
  - "gene prediction rna-seq"
name: braker3-pipeline
description: Run BRAKER3 for eukaryotic gene prediction. Supports RNA-seq + protein, protein-only, and RNA-seq-only modes. Covers alignment, execution, container setup, and BRAKER-specific failure modes.
---

# BRAKER3 Pipeline

**"Run BRAKER3 on my soft-masked genome"** -> Predict protein-coding genes using GeneMark-ETP + AUGUSTUS + TSEBRA. BRAKER3 handles all evidence combinations: RNA-seq + protein, protein-only, and RNA-seq-only.

> **Note:** BRAKER3 subsumes the legacy BRAKER1 (RNA-seq only) and BRAKER2 (protein only) modes — they are now execution modes of the same `braker.pl` script.

## Prerequisites

- **Environment**: Active Pixi environment. Refer to `skills/pixi-skill/SKILL.md`.
- **Upstream evidence** (provide at least one):
  - Soft-masked assembly (`.fasta`). See `./repeat-annotation/SKILL.md`.
  - RNA-seq reads (`.fq.gz`) or sorted BAM from `skills/rna-seq-skills/read-mapping-alignment/SKILL.md` — optional but improves accuracy.
  - OrthoDB protein partition (`.fa`) for the smallest clade containing target species — required for protein modes.
- **Assembly QC**: BUSCO (genome mode), N50, contamination check passed.

## Installation

### Via Pixi

```bash
pixi add braker3 hisat2 samtools
# BRAKER3 bundles GeneMark-ETP, AUGUSTUS, and TSEBRA internally
```

### Via Docker (Recommended — avoids GeneMark key issues)

```bash
docker pull teambraker/braker3:latest
# Run via:
# docker run -v $(pwd):/data teambraker/braker3:latest braker.pl --genome=/data/masked.fa ...
```

## Procedure

### 1. Align RNA-seq (if using RNA-seq and not already aligned)

Use a splice-aware aligner. Output sorted, indexed BAM.

```bash
hisat2-build masked.fasta idx
hisat2 -x idx -1 R1.fq.gz -2 R2.fq.gz --dta -p 16 | samtools sort -@4 -o rnaseq.bam
samtools index rnaseq.bam
```

### 2. Run BRAKER3 — Full Evidence (RNA-seq + Protein)

```bash
braker.pl \
  --genome=masked.fasta \
  --prot_seq=OrthoDB/Vertebrata.fa \
  --bam=rnaseq.bam \
  --softmasking \
  --threads=16 \
  --species=my_species \
  --gff3 \
  --workingdir=braker3_out
```

**Auto-download mode** (BRAKER3 can fetch and align SRA reads itself):

```bash
braker.pl \
  --genome=masked.fasta \
  --prot_seq=OrthoDB/Vertebrata.fa \
  --rnaseq_sets_ids=SRR12345678,SRR12345679 \
  --rnaseq_sets_dirs=/path/to/fastq/ \
  --softmasking \
  --threads=16 \
  --species=my_species \
  --gff3 \
  --workingdir=braker3_out
```

### 2b. Run BRAKER3 — Protein-Only Mode (No RNA-seq)

When RNA-seq is unavailable, run BRAKER3 with **only** `--prot_seq`. This triggers the legacy BRAKER2 protein-homology pipeline internally (GeneMark-EP+ / ProtHint).

```bash
braker.pl \
  --genome=masked.fasta \
  --prot_seq=OrthoDB/Vertebrata.fa \
  --softmasking \
  --threads=16 \
  --species=my_species \
  --gff3 \
  --workingdir=braker3_protein_only_out
```

> **Key point:** The same `braker.pl` script handles protein-only mode — no need for a separate BRAKER2 installation. Simply **omit `--bam` and `--rnaseq_sets_ids`**.

### 2c. Run BRAKER3 — RNA-seq-Only Mode (No Protein)

When protein evidence is unavailable, run with **only** `--bam` or `--rnaseq_sets_ids`. This triggers the legacy BRAKER1 pipeline internally (GeneMark-ET).

```bash
braker.pl \
  --genome=masked.fasta \
  --bam=rnaseq.bam \
  --softmasking \
  --threads=16 \
  --species=my_species \
  --gff3 \
  --workingdir=braker3_rnaseq_only_out
```

### 3. Key Agentic Flags

| Flag | Meaning | Required? |
|------|---------|-----------|
| `--genome` | Soft-masked assembly | **Mandatory** |
| `--prot_seq` | OrthoDB clade partition (smallest containing target) | Required for protein modes |
| `--bam` | Pre-aligned sorted RNA-seq BAM | Required for RNA-seq modes |
| `--rnaseq_sets_ids` | Comma-separated SRA accessions (auto-download) | RNA-seq mode alternative |
| `--softmasking` | Honor lowercase as soft-mask | **Mandatory** |
| `--species` | Unique species name for AUGUSTUS parameter set | **Mandatory** |
| `--gff3` | Emit GFF3 in addition to GTF | Optional |
| `--workingdir` | Output directory | Optional (defaults to current) |

### 4. Outputs

| File | Content |
|------|---------|
| `braker.gtf` / `braker.gff3` | TSEBRA-combined final gene models |
| `braker.codingseq` | Coding sequences (CDS) |
| `braker.aa` | Predicted proteins (FASTA) |

## BRAKER3-Specific Pitfalls

### GeneMark key expiration
**Symptom:** BRAKER fails mid-run with GeneMark license error. **Fix:** Use the official Docker container; recent versions have relaxed this requirement but verify.

### Special characters in FASTA headers
**Symptom:** Parser crashes on spaces or pipes in `>contig_1 | length=1000`. **Fix:** `sed 's/ .*//' genome.fa` to strip descriptions.

### Wrong OrthoDB partition
**Symptom:** Clean-looking GFF3 but systematically mis-placed exons. **Fix:** Pick the smallest OrthoDB partition containing the target clade (e.g., Vertebrata for a fish, not Metazoa).

### Running without `--softmasking`
**Symptom:** 2x inflated gene count, thousands of spurious single-exon genes inside TEs. **Fix:** Always pass `--softmasking` with a soft-masked assembly.

### Expecting isoforms/UTRs
**Symptom:** mRNA:gene ratio == 1.00; downstream AS/APA/scRNA-seq analyses fail. **Fix:** BRAKER3 emits one CDS-only model per locus. Add a PASA/Iso-Seq update step if isoforms are needed.

## Quantitative Checks

| Metric | Expected | Flag if |
|--------|----------|---------|
| Gene count | Within 1.5× of nearest relative (adjust for ploidy) | >2× or <0.5× |
| Mono-exonic fraction | ~10-20% (vertebrate); fungi higher | >30% |
| Protein length | Unimodal ~300-450 aa | Spike <100 aa |
| BUSCO-Duplicated | 1-3% (clean haploid) | >5-8% without known WGD |
| mRNA:gene ratio | >1.00 (if PASA updated) | == 1.00 (isoform-naive) |

## Verification

- [ ] `--softmasking` was passed and genome was soft-masked upstream.
- [ ] Assembly QC passed before prediction (BUSCO genome, N50, no contamination).
- [ ] OrthoDB partition is the smallest clade containing target species.
- [ ] BRAKER completed without GeneMark key or header errors.
- [ ] Output `braker.gff3` sanity-checked (gene count, mono-exonic fraction, mRNA:gene ratio).
- [ ] Protein-BUSCO score vs genome-BUSCO score shows no major gene drop.
- [ ] If isoforms/UTRs needed, PASA/Iso-Seq update step is planned or executed.

## Related Skills

- `./gene-prediction-decision-tree/SKILL.md` — Router: choose the right tool
- `./repeat-annotation/SKILL.md` — PREREQUISITE: soft-mask assembly
- `skills/rna-seq-skills/read-mapping-alignment/SKILL.md` — Splice-aware RNA-seq alignment
- `./annotation-qc/SKILL.md` — Assess annotation quality post-prediction
- `./galba-pipeline/SKILL.md` — Protein-only alternative (close relatives)

## References

- Gabriel L, et al. 2024. BRAKER3: fully automated genome annotation using RNA-seq and protein evidence with GeneMark-ETP, AUGUSTUS, and TSEBRA. *Genome Res* 34:769-777.
- Brůna T, et al. 2024. GeneMark-ETP significantly improves the accuracy of automatic annotation of large eukaryotic genomes. *Genome Res* 34:757-768.
- Gabriel L, et al. 2021. TSEBRA: transcript selector for BRAKER. *BMC Bioinformatics* 22:566.
- Hoff KJ, et al. 2016. BRAKER1: unsupervised RNA-Seq-based genome annotation. *Bioinformatics* 32:767-769.
- Brůna T, et al. 2021. BRAKER2: automatic eukaryotic genome annotation with protein support. *NAR Genom Bioinform* 3:lqaa108.
