---
version: 2
updated: "2026-06-29"
triggers:
  - "gene prediction"
  - "which gene predictor"
  - "gene prediction decision tree"
name: gene-prediction-decision-tree
description: Router skill for eukaryotic gene prediction. Select the optimal tool based on available evidence, assembly quality, and clade. Covers tool taxonomy, decision tree, quantitative thresholds, and cross-cutting failure modes.
---

# Gene Prediction Decision Tree

**"Which gene prediction tool should I use?"** -> Choose the optimal pipeline based on your evidence type, assembly quality, and target clade.

## The Single Most Important Modern Insight

Gene-prediction accuracy is governed by the quality of the **training set** and **extrinsic evidence**, not by which algorithm is chosen. A 2024-era HMM (AUGUSTUS) trained on clean, evidence-validated genes beats a fancier model trained on garbage.

- **The #1 silent killer is training on a bad assembly.** Decontaminate and check assembly contiguity/BUSCO **before** any self-trainer touches the genome.
- BRAKER3's accuracy jump comes from GeneMark-ETP mining a **high-confidence training set** from loci where RNA-seq transcripts AND protein homology independently agree — not from a better HMM.
- The 2025-2026 DL twist: Tiberius/Helixer match BRAKER3 on well-represented clades with **no evidence**, but degrade on under-represented lineages and are isoform/UTR-naive.

## Tool Taxonomy

| Pipeline               | Evidence used                                            | When                                                                             |
| ---------------------- | -------------------------------------------------------- | -------------------------------------------------------------------------------- |
| **BRAKER3**            | RNA-seq + protein, **protein only**, or **RNA-seq only** | **Default for all evidence combinations**; subsumes legacy BRAKER1 and BRAKER2   |
| **GALBA**              | protein only (close relatives)                           | Beats BRAKER3 protein-only on large vertebrate genomes with good close proteomes |
| **GeMoMa**             | reference annotation                                     | Project existing close-relative annotation                                       |
| **Funannotate**        | any                                                      | **Fungal de facto standard**                                                     |
| **Helixer / Tiberius** | none (DL ab initio)                                      | Evidence-poor genomes in well-represented clades                                 |
| **Prokka / Bakta**     | prokaryotic                                              | **Bacterial/archaeal standard**; not covered in this pipeline                    |

OrthoDB rule: BRAKER3 wants a **clade partition** (smallest containing target clade). GALBA wants a few **close-relative proteomes**.

> **Prokaryotic genomes:** This pipeline and all its predictors (BRAKER3, GALBA, Funannotate, Helixer) are designed for **eukaryotic** gene prediction. Prokaryotic genomes need Prokka or Bakta, which handle operons, Shine-Dalgarno motifs, and the very different gene structure of bacteria and archaea. The ncRNA annotation (Phase 3) and QC phases are still applicable to prokaryotes — use `tRNAscan-SE -B` (bacterial) or `-A` (archaeal), and `compleasm -l bacteria`.

> **Legacy note:** BRAKER1 (RNA-seq only) and BRAKER2 (protein only) are now execution modes of BRAKER3. The same `braker.pl` script handles all three evidence combinations — simply omit the evidence flags you don't have.

## Decision Tree by Scenario

| Evidence / scenario                    | Recommended                     | Skill to invoke                               |
| -------------------------------------- | ------------------------------- | --------------------------------------------- |
| RNA-seq + protein DB                   | **BRAKER3** (full evidence)     | `./braker3-pipeline/SKILL.md` §2              |
| RNA-seq only                           | **BRAKER3** (RNA-seq-only mode) | `./braker3-pipeline/SKILL.md` §2c             |
| Protein only, close relatives          | **GALBA**                       | `./galba-pipeline/SKILL.md`                   |
| Protein only, broad/distant            | **BRAKER3** (protein-only mode) | `./braker3-pipeline/SKILL.md` §2b             |
| Prokaryote (bacteria/archaea)          | **Prokka / Bakta**              | Not in this pipeline — see note below         |
| No evidence, represented clade         | Helixer / Tiberius              | Decision-tree only (not yet standalone skill) |
| Reference annotation of close relative | GeMoMa                          | Decision-tree only (reactive skill)           |
| Fungus (any evidence)                  | **Funannotate**                 | `./funannotate-pipeline/SKILL.md`             |
| Need isoforms + UTRs                   | Add PASA / Iso-Seq update       | See specific pipeline skill                   |
| Genome not yet masked                  | -> **Soft-mask first**          | `./repeat-annotation/SKILL.md`                |
| Want to assess the result              | -> **Annotation QC**            | `./annotation-qc/SKILL.md`                    |

> **BRAKER3 mode summary:** Same script, three modes:
> 
> - Full: `--prot_seq` + `--bam`
> - Protein-only: `--prot_seq` only (omit `--bam`)
> - RNA-seq-only: `--bam` only (omit `--prot_seq`)

## Prerequisites

- **Environment**: Active Pixi environment. Refer to `skills/pixi-skill/SKILL.md`.
- **Mandatory prerequisite**: Soft-masked assembly. See `./repeat-annotation/SKILL.md`.
- **Assembly QC**: BUSCO (genome mode), N50, contamination check via FCS-GX/BlobTools.

## Quantitative Thresholds

| Threshold                                                 | Rationale                                                               |
| --------------------------------------------------------- | ----------------------------------------------------------------------- |
| Soft-mask before prediction                               | Unmasked TEs inflate spurious genes                                     |
| Gene count vs nearest relative (±, reconcile with ploidy) | 1.5-2x with no WGD = haplotigs/over-prediction; ~0.5x = over-masking    |
| Mono-exonic fraction ~10-20% (vertebrate)                 | >25-30% = unmasked TEs/pseudogenes/fragments; fungi legitimately higher |
| Protein length unimodal ~300-450 aa                       | Sub-100-aa spike = spurious/fragmented; fat left tail = partials        |
| BUSCO-Duplicated ~1-3% (clean haploid)                    | >5-8% with no WGD -> purge_dups before annotating                       |
| mRNA:gene ratio == 1.00                                   | Isoform/UTR-naive; AS analyses untrustworthy                            |
| Pick smallest OrthoDB partition containing the clade      | Broader = noisier hints, slower                                         |

## Gene-Model Sanity Statistics

**Goal:** Compute the triage panel that reveals annotation health where gene count and BUSCO cannot — isoform ratio, mono-exonic fraction, protein-length distribution.

```python
import gffutils

MONOEXONIC_FLAG = 0.30   # >30% single-exon in a vertebrate suggests unmasked TEs/pseudogenes/fragments

def gene_model_stats(gff_file):
    db = gffutils.create_db(gff_file, ':memory:', merge_strategy='merge')
    genes = list(db.features_of_type('gene'))
    mrnas = list(db.features_of_type(['mRNA', 'transcript']))
    exon_counts = [len(list(db.children(tx, featuretype='exon'))) for tx in mrnas]
    mono_frac = sum(1 for e in exon_counts if e == 1) / len(exon_counts) if exon_counts else 0
    mrna_per_gene = len(mrnas) / len(genes) if genes else 0
    if mrna_per_gene <= 1.001:
        print('WARNING: one isoform per locus (mRNA:gene == 1.00) -- isoform/UTR-naive')
    if mono_frac > MONOEXONIC_FLAG:
        print(f'WARNING: mono-exonic fraction {mono_frac:.1%} high -- check masking/contamination')
    return {'genes': len(genes), 'mrna_per_gene': mrna_per_gene, 'mono_exonic_fraction': mono_frac}
```

## Cross-Cutting Pitfalls

### Unmasked or hard-masked genome

**Symptom:** 2x inflated gene count, high mono-exonic fraction, fragmented models. **Fix:** soft-mask with curated library; pass `--softmasking`.

### Training on a bad assembly

**Symptom:** Confidently wrong, BUSCO-green models genome-wide. **Fix:** decontaminate and check N50/BUSCO first.

### Wrong OrthoDB partition

**Symptom:** Systematically mis-placed exon boundaries; clean-looking GFF3. **Fix:** pick smallest clade partition containing target species.

### Expecting isoforms/UTRs from a one-model pipeline

**Symptom:** Discarded scRNA-seq reads; empty AS results. **Fix:** add PASA/Iso-Seq update step; check mRNA:gene ratio.

### High BUSCO-Duplicated read as success

**Symptom:** Inflated gene count. **Fix:** if no known WGD and D>5-8%, purge_dups first.

## Common Errors

| Error / symptom                      | Cause                                              | Solution                                 |
| ------------------------------------ | -------------------------------------------------- | ---------------------------------------- |
| Thousands of extra short genes       | Unmasked repeats                                   | Soft-mask; `--softmasking`               |
| BRAKER fails mid-run                 | Expired GeneMark key / special chars in headers    | Use container; `sed 's/ .*//' genome.fa` |
| Many single-exon genes               | Unmasked TEs / contamination / fragmented assembly | Verify masking; decontaminate; check N50 |
| Low protein-BUSCO, high genome-BUSCO | Predictor missed present genes                     | Fix evidence/masking, not assembly       |
| AS analysis returns nothing          | One-isoform annotation                             | Run PASA/Iso-Seq update                  |
| Suspiciously few NLR/ZNF genes       | Over-masked with uncurated library                 | Filter repeat library against protein DB |

## Verification

- [ ] Genome is soft-masked before prediction.
- [ ] Assembly QC passed (BUSCO, N50, contamination check).
- [ ] Correct tool selected by evidence type.
- [ ] OrthoDB partition is the smallest clade containing target species.
- [ ] Output GFF3 sanity-checked (gene count, mono-exonic fraction, mRNA:gene ratio).
- [ ] If isoform/UTR analysis needed, PASA/Iso-Seq update step planned.
- [ ] Protein-BUSCO vs genome-BUSCO comparison shows no major missing gene issue.

## Related Skills

- `./repeat-annotation/SKILL.md` — PREREQUISITE: soft-mask repeats
- `./braker3-pipeline/SKILL.md` — Default for all evidence combinations (RNA-seq + protein, protein-only, RNA-seq-only)
- `./galba-pipeline/SKILL.md` — Protein-only close relatives
- `./funannotate-pipeline/SKILL.md` — Fungi
- `./annotation-qc/SKILL.md` — Assess annotation quality
- `./ncrna-annotation/SKILL.md` — Non-coding RNAs (not found by protein predictors)

## References

- Brůna T, et al. 2024. GeneMark-ETP significantly improves the accuracy of automatic annotation of large eukaryotic genomes. *Genome Res* 34:757-768.
- Gabriel L, et al. 2024. BRAKER3: fully automated genome annotation using RNA-seq and protein evidence. *Genome Res* 34:769-777.
- Brůna T, et al. 2023. GALBA: genome annotation with miniprot and AUGUSTUS. *BMC Bioinformatics* 24:327.
- Palmer JM, Stajich J. 2020. Funannotate v1.8. Zenodo.
- Gabriel L, et al. 2024. Tiberius: end-to-end deep learning with an HMM for gene prediction. *Bioinformatics* 40:btae685.
- Weisman CM, et al. 2020. Many, but not all, lineage-specific genes can be explained by homology detection failure. *PLoS Biol* 18:e3000862.
