---
version: 2
updated: "2026-06-29"
triggers:
  - "repeat annotation"
  - "soft-mask repeats"
  - "earlgrey"
  - "repeatmodeler"
name: repeat-annotation
description: Discover, classify, and soft-mask repetitive elements and transposable elements. Covers RepeatModeler2 (de novo library), RepeatMasker (masking), EDTA (plants), EarlGrey (non-model wrapper), and TE expression from RNA-seq. Critical prerequisite for eukaryotic gene prediction.
---

# Repeat and Transposable Element Annotation

**"Mask repeats in my genome assembly"** -> Build a de novo repeat-family library, annotate copies genome-wide, and **soft-mask** them as a mandatory prerequisite for gene prediction.

## The Single Most Important Modern Insight

Two load-bearing truths the masker hides:

1. **De novo library construction is a curation research project, not a button.** A RepeatModeler2 run emits `mydb-families.fa` overnight — a *draft of a draft*: consensi are routinely 5'-truncated, chimeric, and 30–60% "Unknown" on a non-model genome. Masking percentage is robust to a bad library, but downstream biology rots: inflated family counts, wrong classification, and host-gene-contaminated consensi that silently mask real genes.
2. **Annotation quality is capped by assembly quality.** Short-read assemblies *collapse* near-identical TE copies and *drop* the youngest ones, systematically under-counting TEs. Always ask what assembly a "% repeat" came from.

## Tool Taxonomy

| Tool               | Role                                                        | When                                         |
| ------------------ | ----------------------------------------------------------- | -------------------------------------------- |
| **RepeatModeler2** | de novo family discovery → consensus library                | Discover genome-specific families            |
| **RepeatMasker**   | annotate/mask genome **against** a library                  | The masking step; does not discover families |
| **EarlGrey**       | wraps RepeatModeler2 + auto-curation + RepeatMasker + plots | **Non-model default**; minimal hand-work     |
| **EDTA**           | structural LTR/TIR/Helitron discovery + filtering           | **Plant / structurally-rich** genomes        |
| **LTR_retriever**  | isolate intact LTR-RTs; feeds LAI                           | LTR focus / assembly-quality check           |
| **TRF**            | tandem/satellite repeats                                    | Different algorithm class from TE maskers    |

> **Database schism:** Dfam (open) vs RepBase (paywalled since 2019). Modern open pipelines build on Dfam + de novo. Record the Dfam release and library provenance.

## Decision Tree by Scenario

| Scenario                                   | Recommended                                                | Why                                            |
| ------------------------------------------ | ---------------------------------------------------------- | ---------------------------------------------- |
| Non-model eukaryote, minimal hand-work     | **EarlGrey**                                               | RepeatModeler2 + auto-curation + clean outputs |
| Plant / structurally-rich TE genome        | **EDTA**                                                   | Best-in-class LTR/TIR/Helitron annotation      |
| Well-covered vertebrate, just need masking | **RepeatMasker** `-species` against Dfam                   | Curated families already exist                 |
| Mask before gene prediction                | **RepeatMasker** `-xsmall` (soft) + decontaminated library | Predictors need soft-masking                   |
| Publication-grade TE biology claim         | de novo -> manual curation (TE-Aid, MCHelper)              | Automated library is the start, not the end    |

## Prerequisites

- **Environment**: Active Pixi environment. Refer to `skills/pixi-skill/SKILL.md`.
- **Input**: Assembled genome (`.fa`). Assembly quality caps TE annotation (short-read = systematic undercount).

## Installation

```bash
pixi add repeatmodeler repeatmasker trf
# Optional wrappers:
pixi add edta earlgrey
# For TE expression analysis:
pixi add tetranscripts
```

## Procedure

### 1. De Novo Family Discovery (RepeatModeler2)

```bash
# Build the database
BuildDatabase -name mydb assembly.fa

# Run discovery (-LTRStruct enables LTR structural pipeline)
RepeatModeler -database mydb -threads 16 -LTRStruct
```

Output: `mydb-families.fa` (draft library — expect 30–60% "Unknown" on non-model genomes).

### 2. (Recommended) Decontaminate the Library

BLAST the library against a protein database. Drop any consensus whose best hit is a host gene with no transposase/RT/integrase domain. This prevents the **domesticated-gene massacre** (RAG1, CENP-B, KRAB-ZNF arrays, etc.).

```bash
# Example: filter against UniProt/SwissProt
blastp -query mydb-families.fa -db swissprot -outfmt 6 -evalue 1e-5 > library_vs_proteins.blastp
# Inspect hits; remove gene-family contaminants from the .fa before masking
```

### 3. Union with Dfam Clade (Recommended)

Combine de novo families with the Dfam partition for your clade to maximize sensitivity:

```bash
cat mydb-families.fa dfam_clade_partition.fa > combined_library.fa
```

### 4. Soft-Mask with RepeatMasker

```bash
RepeatMasker -lib combined_library.fa -xsmall -gff -e rmblast -pa 16 -dir rm_out assembly.fa
```

**Key flags:**
| Flag | Meaning |
|------|---------|
| `-lib` | Custom library (de novo + Dfam union) |
| `-xsmall` | **Soft-mask (lowercase)** — mandatory for gene prediction |
| `-gff` | Emit GFF annotation |
| `-e rmblast` | Engine: consensus FASTA alignment (default) |
| `-pa 16` | Parallel threads |
| `-nolow` | Skip low-complexity/simple repeats (recommended before gene prediction) |

> **Critical:** Default RepeatMasker hard-masks with `N`. `-xsmall` gives lowercase. For gene prediction, also use `-nolow` because simple repeats overlap real coding microsatellites.

### 5. Outputs

| File                 | Content                                  |
| -------------------- | ---------------------------------------- |
| `assembly.fa.masked` | Soft-masked genome (lowercase = repeat)  |
| `assembly.fa.out`    | Tabular annotation of each repeat copy   |
| `assembly.fa.tbl`    | Summary statistics (% by class)          |
| `assembly.fa.align`  | Alignments (needed for Kimura landscape) |

## Soft vs Hard Masking (Critical)

| Masking type  | Flag      | Result              | When to use                                                                               |
| ------------- | --------- | ------------------- | ----------------------------------------------------------------------------------------- |
| **Soft-mask** | `-xsmall` | Repeats → lowercase | **Gene prediction** — predictor avoids nucleating in repeats but exons can extend through |
| Hard-mask     | default   | Repeats → `N`       | **NEVER for gene prediction** — truncates genes overlapping repeats silently              |
| X-mask        | `-x`      | Repeats → `X`       | Rarely used                                                                               |

**The over-masking massacre:** An uncurated de novo library contains fragments of real multicopy gene families. Masking them deletes the genome's most interesting genes: RAG1/RAG2, syncytins, SETMAR/Metnase, CENP-B, KRAB-ZNF arrays. **Decontaminate the library against a protein DB before masking.**

## Repeat Statistics and Age Landscape

Parse the `.out` file to summarize by class and plot Kimura divergence (~relative age).

```python
import pandas as pd
import matplotlib.pyplot as plt

def parse_repeatmasker_out(out_file):
    records = []
    with open(out_file) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            records.append({
                'perc_div': float(parts[1]),
                'seqid': parts[4],
                'repeat_class': parts[10],
                'length': int(parts[6]) - int(parts[5]) + 1
            })
    return pd.DataFrame(records)

def repeat_summary(rm_df, genome_size):
    by_class = rm_df.groupby('repeat_class')['length'].sum().sort_values(ascending=False)
    total = rm_df['length'].sum()
    print(f'Total masked: {total/genome_size:.1%} of genome')
    print('Note: this is a LOWER bound; ancient copies decay past detection')
    return by_class / genome_size * 100

def kimura_landscape(rm_df, bins=50):
    """Plot divergence histogram by class. Low divergence = recent burst."""
    classes = rm_df['repeat_class'].unique()
    fig, ax = plt.subplots()
    for cls in classes[:6]:  # Top 6 classes
        subset = rm_df[rm_df['repeat_class'] == cls]['perc_div']
        ax.hist(subset, bins=bins, alpha=0.6, label=cls)
    ax.set_xlabel('Kimura divergence from consensus (%)')
    ax.set_ylabel('Copy count')
    ax.legend()
    return fig
```

> **Caution:** The landscape is **right-censored** — ancient TEs decayed past detection, so "no old activity" can mean "invisible," not "absent."

## TE Expression from RNA-seq (Multimapping Minefield)

Young high-copy TE families produce reads that map equally to hundreds of near-identical loci. **Unique-only mapping** discards most TE signal and biases toward old, uniquely-mappable copies.

**Solution:** EM/probabilistic reassignment tools:

- **TEtranscripts / TElocal** — subfamily-level (high power, no locus resolution)
- **SQuIRE** — locus-level (noisy, mappability-sensitive)
- **Telescope** — locus-level probabilistic

**Dominant false positive:** A TE in an intron or downstream of an expressed gene is **not** "expressed" — it's read-through. Distinguish autonomous transcription from passenger signal by strand and continuity.

## Pitfalls

### Hard-masking before gene prediction

**Symptom:** Genes overlapping repeats silently absent from GFF. **Fix:** Always use `-xsmall` for gene prediction.

### Host-gene-contaminated library

**Symptom:** Suspiciously few NLR/ZNF/OR/RAG1 genes. **Fix:** Decontaminate library against protein DB before masking.

### Trusting % repeat from a short-read assembly

**Symptom:** "Low TE, all ancient" or non-comparable cross-study tables. **Fix:** Check LAI/assembly type; report method + assembly with every number.

### Discarding multimappers in TE RNA-seq

**Symptom:** Most TE signal lost, bias to old elements. **Fix:** Use EM tools (TEtranscripts/SQuIRE/Telescope).

### "Unknown" passed off as a result

**Symptom:** 40% Unknown library shipped without inspection. **Fix:** Run RepeatClassifier/DeepTE to triage; curate; note DB-coverage limits.

## Quantitative Thresholds

| Threshold                                      | Rationale                                                                         |
| ---------------------------------------------- | --------------------------------------------------------------------------------- |
| `-xsmall` soft-mask before gene prediction     | Hard-mask truncates repeat-overlapping genes                                      |
| TE content scales with genome size             | Human ~50%, maize ~85%, Arabidopsis ~20–25%, fungi ~1–20% — sanity-check vs clade |
| "Unknown" ~<15% (mammal) vs 30–50% (non-model) | High Unknown bounds biological claims                                             |
| LAI <10 draft / 10–20 reference / >20 gold     | LTR-RT resolution metric; only valid for LTR-rich genomes                         |
| Report library + engine + assembly with any %  | % masked is non-comparable across methods                                         |

## Verification

- [ ] Library decontaminated against protein DB (no domesticated-gene massacre).
- [ ] `-xsmall` used for soft-masking (not default hard-mask).
- [ ] `-nolow` considered before gene prediction (preserves coding microsatellites).
- [ ] RepeatMasker completed without errors.
- [ ] `assembly.fa.masked` exists and contains lowercase letters.
- [ ] `assembly.fa.tbl` shows expected % for clade/genome size.
- [ ] If publication-grade claims: library manually curated or EarlGrey/EDTA used.
- [ ] If TE expression planned: EM multimapping tool selected (not unique-only).

## Related Skills

- `./gene-prediction-decision-tree/SKILL.md` — Receives the soft-masked genome
- `./braker3-pipeline/SKILL.md` — Gene prediction requires soft-masked input
- `./annotation-qc/SKILL.md` — LAI and repeat-content sanity checks
- `skills/rna-seq-skills/differential-expression-analysis/` — Differential TE expression from TEtranscripts counts (sub-skill of the rna-seq-skills meta-skill)

## References

- Flynn JM, et al. 2020. RepeatModeler2 for automated genomic discovery of transposable element families. *PNAS* 117:9451–9457.
- Storer J, et al. 2021. The Dfam community resource of transposable element families. *Mob DNA* 12:2.
- Ou S, et al. 2019. Benchmarking transposable element annotation methods (EDTA). *Genome Biol* 20:275.
- Ou S, Jiang N. 2018. LTR_retriever: identification of long terminal repeat retrotransposons. *Plant Physiol* 176:1410–1422.
- Baril T, et al. 2024. Earl Grey: a fully automated transposable element annotation pipeline. *Mol Biol Evol* 41:msae068.
- Goubert C, et al. 2022. A beginner's guide to manual curation of transposable elements. *Mob DNA* 13:7.
- Jin Y, et al. 2015. TEtranscripts: including transposable elements in DE analysis. *Bioinformatics* 31:3593–3599.
- Wicker T, et al. 2007. A unified classification system for eukaryotic transposable elements. *Nat Rev Genet* 8:973–982.
