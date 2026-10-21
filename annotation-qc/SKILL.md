---
version: 2
updated: "2026-06-29"
triggers:
  - "annotation QC"
  - "BUSCO proteome"
  - "OMArk"
  - "is my annotation any good"
name: annotation-qc
description: Assess genome annotation quality with BUSCO, compleasm (fast miniprot-based alternative), OMArk, CheckM2, and gene-set sanity metrics. Covers the assembly-vs-proteome diagnostic, what duplicated scores mean, and why gene count is not a quality metric.
---

# Annotation Quality Control

**"Is my genome annotation any good?"** -> Measure conserved-gene completeness, proteome consistency/contamination, and gene-set sanity, and decide whether the limiting factor is the assembly or the predictor.

## The Single Most Important Modern Insight

BUSCO completeness measures recovery of the ~1,000–5,500 most conserved single-copy orthologs — the genes any half-broken pipeline will find. A 98% BUSCO confirms the housekeeping core is intact; it says **nothing** about whether the other ~20,000 models are chimeric, fragmented, fused, or hallucinated from TEs. BUSCO is a smoke detector in one room of a burning house.

**Three load-bearing moves:**
1. **Run assembly-BUSCO vs proteome-BUSCO on the same assembly — the diagnostic almost nobody runs.** `-m genome` does its own gene-finding; `-m proteins` scores the delivered proteome. If assembly is 98% and proteome is 85%, the genes are present and the *predictor* missed them — fix training/evidence/masking. If both are low, stop annotating and fix the assembly.
2. **Gene count is a vanity metric.** The same assembly annotated by two pipelines routinely differs 20–40% in gene count. The diagnostic signal is in the **mono-exonic fraction, protein-length distribution, mean exons/gene, and mRNA:gene ratio** — not the headline count.
3. **Complement BUSCO with OMArk.** BUSCO's single-copy lens is blind to over-prediction, chimeras, and contamination; OMArk assesses the *whole* proteome for completeness AND consistency.

## Tool Taxonomy

| Tool | Measures | Domain | When |
| --- | --- | --- | --- |
| **BUSCO** | Conserved single-copy ortholog recovery (C/D/F/M) | Euk + prok + viral | Field-standard; conservative |
| **compleasm** | Faster miniprot-based BUSCO-compatible completeness | Euk + prok | Speed-critical; large genomes; often higher completeness than metaeuk |
| **OMArk** | Proteome completeness + consistency + contamination | Eukaryotic proteomes | Catches over-prediction BUSCO misses |
| **CheckM2** | Completeness + contamination (ML, lineage-agnostic) | Bacteria/archaea | Prokaryotic standard; handles MAGs |
| **Gene-set sanity** | Gene count, mono-exonic %, protein length, mRNA:gene, coding density | Any | The panel BUSCO can't provide |
| **LAI** | LTR-RT assembly resolution | LTR-rich genomes | Assembly-side QC before annotation |

## Decision Tree by Scenario

| Scenario | QC to run | Why |
| --- | --- | --- |
| Eukaryotic annotation (standard) | BUSCO `-m proteins` + compleasm `run` + OMArk + sanity panel | Genome-vs-proteome fork + over-prediction check |
| Eukaryotic annotation (speed-critical / large genome) | **compleasm** `run` + OMArk + sanity panel | compleasm is ~10–100× faster than BUSCO; same lineage DBs |
| Prokaryotic genome/MAG | CheckM2 + coding density + tRNA/rRNA counts | CheckM2 is the field-standard for prokaryotes |
| Suspect high gene count | Mono-exonic fraction + protein-length + BUSCO/compleasm-D | Distinguish haplotigs / unmasked TEs / split models |
| High duplicated score | Check synteny + clade ploidy | Uncollapsed haplotigs (purge) vs real WGD (keep) |
| Low proteome-BUSCO, high genome-BUSCO | Fix evidence/masking/training | Genes are present; predictor missed them |
| Both BUSCO/compleasm low | -> Assembly QC (purge_dups, more data) | Limit is the assembly |

## Prerequisites

- **Environment**: Active Pixi environment. Refer to `skills/pixi-skill/SKILL.md`.
- **Input**: Predicted proteins (`.faa`) and/or genome assembly (`.fa`), annotation GFF3.
- **Databases**: BUSCO/compleasm lineage dataset (e.g. `vertebrata_odb10` or `odb12`), OMArk LUCA database.

## Installation

```bash
pixi add busco omark checkm2 compleasm gffutils pandas matplotlib
```

Download OMArk database (do once):
```bash
omamer download --db LUCA.h5
```

## Procedure

### 1. BUSCO — Field-Standard Completeness

```bash
# Proteome mode — the delivered annotation
busco -i proteins.faa -m proteins -l vertebrata_odb10 -o busco_prot -c 16

# Genome mode — assembly capacity (for the diagnostic comparison)
busco -i genome.fa -m genome -l vertebrata_odb10 -o busco_genome -c 16
```

**Output categories:** `C:[S,D],F,M` — Complete (Single/Duplicated), Fragmented, Missing.

**Use the deepest applicable clade dataset**, not shallow `eukaryota_odb10` (~255 genes). A 99% on 255 genes is trivially easy and not comparable to 99% on 5,500 genes.

### 2. compleasm — Fast Miniprot-Based Alternative

```bash
# Genome mode (the most common use — directly from assembly)
compleasm run -a genome.fa -o compleasm_out -l vertebrata -t 16

# Protein mode (assess delivered proteome, analogous to BUSCO -m proteins)
compleasm protein -p proteins.faa -o compleasm_prot_out -l vertebrata -t 16

# Auto-detect lineage (no need to guess the clade)
compleasm run -a genome.fa -o compleasm_auto -t 16 --autolineage
```

**Output categories:** `S` (Single), `D` (Duplicated), `F` (Fragmented, subclass 1), `I` (Fragmented, subclass 2 — split across contigs), `M` (Missing).

**Key differences from BUSCO:**

| Feature | BUSCO | compleasm |
|---------|-------|-----------|
| Engine | metaeuk (HMM + ORF) | miniprot (protein-to-genome alignment) |
| Speed | Slow (hours for large genomes) | **Fast** (minutes) |
| Fragmented | Single `F` category | `F` + `I` (split across contigs) |
| Completeness | Often lower (conservative) | Often **higher** (miniprot more sensitive) |
| Retrocopy handling | Duplicated (`D`) | `--retrocopy` separates retrocopies from true duplicates |
| Lineage DB | `odb10` / `odb12` | Same `odb10` / `odb12` (default: `odb12` in v0.2.7+) |

**When to prefer compleasm:**
- Large genomes where BUSCO runtime is prohibitive
- Assembly QC pipelines where speed matters
- When you want the assembly-vs-predictor diagnostic quickly
- The `--retrocopy` flag is useful for distinguishing WGD from retrocopies

**When to prefer BUSCO:**
- Conservative, field-standard benchmark for publications
- When reviewers specifically request "BUSCO" by name
- `metaek` path may be more sensitive for highly divergent organisms

### 3. OMArk — Proteome Consistency and Contamination

```bash
omamer search --db LUCA.h5 --query proteins.faa --out proteins.omamer
omark -f proteins.omamer -d LUCA.h5 -o omark_out
```

OMArk classifies the *whole* proteome as consistent / inconsistent (structure conflicts with homologs = chimera or fragment) / unknown, and detects **contaminant species**. Use alongside BUSCO/compleasm — it catches over-prediction the single-copy lens misses.

### 4. CheckM2 — Prokaryotic Standard

```bash
checkm2 predict --input genome.fna --output-directory checkm2_out --threads 16
```

**Gate annotation on it:** contamination >5% mixes organisms; completeness <90% with contamination >5–10% makes gene count uninterpretable — fix assembly/binning first.

### 5. Gene-Set Sanity Panel

Load the GFF3 into gffutils and compute the triage metrics BUSCO cannot provide:

```python
import gffutils

MONOEXONIC_FLAG = 0.30   # >30% single-exon in vertebrate = unmasked TEs/pseudogenes

def sanity_panel(gff_file):
    db = gffutils.create_db(gff_file, ':memory:', merge_strategy='merge')
    genes = list(db.features_of_type('gene'))
    mrnas = list(db.features_of_type(['mRNA', 'transcript']))
    exon_counts = [len(list(db.children(tx, featuretype='exon'))) for tx in mrnas]
    mono_frac = sum(1 for e in exon_counts if e == 1) / len(exon_counts) if exon_counts else 0
    mrna_per_gene = len(mrnas) / len(genes) if genes else 0
    if mrna_per_gene <= 1.001:
        print('WARNING: mRNA:gene == 1.00 -- isoform/UTR-naive')
    if mono_frac > MONOEXONIC_FLAG:
        print(f'WARNING: mono-exonic fraction {mono_frac:.1%} high -- check masking/contamination')
    return {'genes': len(genes), 'mrna_per_gene': mrna_per_gene,
            'mono_exonic_fraction': mono_frac}
```

Read gene count against the **nearest well-annotated relative** and ploidy, never in isolation. A healthy protein-length distribution is unimodal near ~300–450 aa; a sub-100-aa spike = spurious/fragmented calls; a fat left tail = partials from a fragmented assembly.

## Pitfalls

### Reporting genome-mode BUSCO as the annotation's score
**Symptom:** Cited BUSCO higher than real proteome quality. **Fix:** Always run `-m proteins` (or `compleasm protein`) on the delivered proteome; report both and compare.

### High Duplicated read as success
**Symptom:** Inflated gene count. **Three causes, opposite responses:**
- **Uncollapsed haplotigs** (genome-wide, non-syntenic) → `purge_dups` *before* annotating
- **Real WGD** (syntenic, Ks peak) → keep
- **Split models** from fragmented assembly → fix assembly

### Shallow lineage dataset
**Symptom:** `eukaryota_odb10` (~255 genes) gives misleadingly high 99%. **Fix:** Use deepest applicable clade; record the dataset name.

### BUSCO alone, no OMArk
**Symptom:** Green BUSCO on a proteome full of chimeras or contamination. **Fix:** Add OMArk for consistency + contamination check.

### Prokaryote: annotating before CheckM2
**Symptom:** Chimeric gene set, uninterpretable metrics. **Fix:** CheckM2 first; contamination >5% → decontaminate.

## Quantitative Thresholds

| Threshold | Rationale |
| --- | --- |
| Assembly-BUSCO vs proteome-BUSCO (or compleasm) gap | Large gap = predictor missed present genes; both low = fix assembly |
| BUSCO-Duplicated ~1–3% (clean haploid); >5–8% no WGD | `purge_dups` before annotating |
| Use deepest applicable clade `_odb10` / `_odb12` | Shallow sets trivially hit 99% |
| Mono-exonic ~10–20% (vertebrate); >25–30% flag | Unmasked TEs/pseudogenes/fragments |
| Protein length unimodal ~300–450 aa | Sub-100-aa spike = spurious; fat left tail = partials |
| mRNA:gene ratio == 1.00 | Isoform/UTR-naive |
| CheckM2 contamination ≤5%, completeness ≥90% | Above/below = prokaryotic QC numbers uninterpretable |
| Prokaryotic coding density ~88–90% | <85% = wrong table/fragmentation; >93% = over-calling |

## Verification

- [ ] **Proteome-mode** BUSCO or compleasm run on delivered proteins (not genome mode alone).
- [ ] Assembly-mode also run for the diagnostic comparison.
- [ ] Deepest applicable clade dataset used (record name: `vertebrata_odb10`, `odb12`, etc.).
- [ ] OMArk run for consistency + contamination check (eukaryotes).
- [ ] CheckM2 run for prokaryotes (completeness >90%, contamination <5%).
- [ ] Gene-set sanity panel computed (mono-exonic %, protein length, mRNA:gene).
- [ ] Duplicated score interpreted against ploidy/synteny (not celebrated blindly).
- [ ] If compleasm used: version and OrthoDB version (`--odb`) recorded.

## Related Skills

- `./braker3-pipeline/SKILL.md` — The annotation this QC evaluates
- `./galba-pipeline/SKILL.md` — Alternative gene predictor to compare
- `./repeat-annotation/SKILL.md` — LAI and over-masking diagnosis
- `./ncrna-annotation/SKILL.md` — tRNA/rRNA count sanity
- `./gene-prediction-decision-tree/SKILL.md` — Router for choosing predictors

## References

- Manni M, et al. 2021. BUSCO update. *Mol Biol Evol* 38:4647–4654.
- Simão FA, et al. 2015. BUSCO: assessing genome assembly and annotation completeness. *Bioinformatics* 31:3210–3212.
- Huang N, Li H. 2023. compleasm: a faster and more accurate reimplementation of BUSCO. *Bioinformatics* 39:btad595.
- Nevers Y, et al. 2025. Quality assessment of gene repertoire annotations with OMArk. *Nat Biotechnol* 43:124–133.
- Chklovski A, et al. 2023. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality. *Nat Methods* 20:1203–1212.
- Guan D, et al. 2020. Identifying and removing haplotypic duplication in primary genome assemblies (purge_dups). *Bioinformatics* 36:2896–2898.
- Ou S, Chen J, Jiang N. 2018. Assessing genome assembly quality using the LTR Assembly Index (LAI). *Nucleic Acids Res* 46:e126.
