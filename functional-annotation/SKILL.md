---
version: 2
updated: "2026-06-29"
triggers:
  - "functional annotation"
  - "GO KEGG Pfam annotation"
  - "eggnog-mapper"
  - "interproscan"
name: functional-annotation
description: Assign GO terms, Pfam/InterPro domains, KEGG orthologs, EC numbers, and product names to predicted proteins using eggNOG-mapper, InterProScan, and KofamScan. Covers annotation-error propagation, domain-presence-is-not-function, GO IEA circularity, and evidence tiering. Use when functionally annotating predicted genes or interpreting gene-family evolution results.
---

# Functional Annotation

**"What do these predicted proteins do?"** -> Transfer GO/KEGG/Pfam/EC/product labels from characterized proteins by orthology and domain signatures, attaching a confidence tier and provenance to each.

## The Single Most Important Modern Insight — Annotation Is a Propagated Hypothesis

Almost every label on a new genome is *transferred* by homology from a small island of experimentally characterized proteins. The transfer chain is lossy and self-reinforcing — a **percolation cascade** (Gilks 2002). An over-specific name assigned once becomes the nearest hit for the next genome, compounding error.

**Three consequences:**
1. **Goal = "honestly tiered," not "maximally annotated."** 40% "hypothetical protein" with correct tiering is better than 95% named with half wrong.
2. **Domain presence ≠ function known.** A Pfam hit reports architecture, not activity. ~10% of human kinome are catalytically dead pseudokinases carrying a confident "protein kinase" domain.
3. **Record provenance on every label.** Method, donor, evidence code, identity/coverage, DB version — the only thing that stops provenance amnesia.

## Tool Taxonomy

| Paradigm | Tool | Mechanism | Failure mode |
|----------|------|-----------|--------------|
| Orthology | **eggNOG-mapper** | Seed ortholog → OG → consensus transfer | Tax-scope sensitive; HGT/xenologs break assumption |
| Domain | **InterProScan** | Profile HMMs → integrated InterPro entries | Domain implies capability, not substrate |
| KEGG ortholog | **KofamScan** | Per-KO HMMs + adaptive thresholds | KO assignment, not pathway proof |
| Best-hit | DIAMOND vs Swiss-Prot | Top-hit similarity | Transitive error propagation |
| ML / structure | DeepGO, Foldseek | Learned sequence/structure → GO | Low precision; reaches twilight zone only |

**Default pair:** eggNOG-mapper + InterProScan (orthogonal evidence: orthology vs signatures). Add KofamScan if KEGG pathway reconstruction is the goal.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Bacterial isolate | Bakta/PGAP + eggNOG + InterProScan | Structural pipeline first, then orthology + domains |
| Eukaryotic proteome | InterProScan + eggNOG-mapper | Orthogonal evidence, reconcile |
| Metagenome / MAG | eggNOG `--itype metagenome` + KofamScan | Built-in gene calling; KEGG modules |
| Twilight-zone / ORFan | DeepGOPlus or Foldseek → DeepFRI | Only handle homology-free fraction; low-confidence |
| CAZymes / BGCs / AMR | → dbCAN / antiSMASH / AMRFinderPlus | Generic Pfam gives no substrate/phenotype |
| GO enrichment downstream | → `skills/gene-family-evolution/SKILL.md` §7 | Mind IEA circularity (see Pitfalls below) |

## Prerequisites

- **Environment**: Active Pixi environment. Refer to `skills/pixi-env-mgmt/SKILL.md`.
- **Input**: Protein FASTA from `./braker3-pipeline/SKILL.md` or `./funannotate-pipeline/SKILL.md`.
- **Compute**: eggNOG DB ~44 GB; InterProScan Java bundle tens of GB.

## Installation

```bash
pixi add eggnog-mapper interproscan kofamscan diamond mmseqs2
# eggNOG data: download_eggnog_data.py --data_dir db/ -y
# InterProScan: download member DB bundle separately
```

## Procedure

### 1. eggNOG-mapper — Orthology-Based Annotation

Three stages: (1) **seed-ortholog search** (DIAMOND best-hit — *not* the annotation); (2) **orthology assignment** retrieves fine-grained orthologs; (3) **functional transfer** pools terms across the ortholog set (damps single-entry misannotation).

```bash
download_eggnog_data.py --data_dir db/ -y   # ~44 GB; pin v2.1.15

emapper.py -i proteins.faa --itype proteins -m diamond \
    --tax_scope auto --data_dir db/ --cpu 16 \
    -o annot --output_dir out/
```

**Key parameter:** `--tax_scope` is the single most consequential setting. `auto` lets each seed take its most-informative phylogenetic ceiling. Too broad = distant orthologs, over-generalized function.

**Output columns:** `seed_ortholog`, `eggNOG_OGs`, `COG_category`, `Description`, `Preferred_name`, `GOs`, `EC`, `KEGG_ko`, `PFAMs` (`-` = empty).

### 2. InterProScan — Domain Signature Annotation

```bash
interproscan.sh -i proteins.faa -f TSV,GFF3 -goterms -pa -cpu 16
```

Runs member DBs (Pfam, PANTHER, CDD, SMART, Gene3D, ...) and **integrates overlapping signatures into InterPro entries** (stable IPRxxxxxx). **Report at InterPro-entry level** — the consensus that survives one member DB being wrong.

| Flag | Meaning |
|------|---------|
| `-goterms` | Add interpro2go mapping (IEA/electronic) |
| `-pa` | Reactome/MetaCyc pathways |
| `-dp` | Force local compute (novel/confidential; reproducibility) |

Default uses EBI precalculated lookup (fast, MD5-keyed). For millions of proteins, chunk FASTA into array jobs.

### 3. KofamScan — KEGG Ortholog Assignment

```bash
exec_annotation -f detail-tsv -o kofam_out.tsv \
    --cpu 16 proteins.faa
```

Adaptive per-KO thresholds (Aramaki 2020). An `*` marks above-threshold hits — stricter than eggNOG's bulk `KEGG_ko`.

### 4. Reconcile Multi-Tool Output

**Goal:** Merge eggNOG + InterProScan per protein while preserving provenance. Prefer orthology `Preferred_name` / `Description` for product name; union GO with source tags.

```python
import pandas as pd

def parse_eggnog(path):
    df = pd.read_csv(path, sep='\t', comment='#', header=None)
    cols = ['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs',
            'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name',
            'GOs', 'EC', 'KEGG_ko']
    df.columns = (cols + [f'c{i}' for i in range(len(df.columns)-len(cols))])[:len(df.columns)]
    return df

def best_product_name(row):
    name = row.get('Preferred_name', '-')
    return name if name not in ('-', '', None) else 'hypothetical protein'
```

Use AGAT to graft results onto GFF3. For GO deliverables use GAF (carries evidence code); keep each tool in its own `Dbxref` namespace.

## Pitfalls

### Annotation-error percolation
**Symptom:** Confident specific names with no provenance. **Fix:** orthology consensus + Swiss-Prot donors; never use DIAMOND-vs-TrEMBL alone.

### Domain-presence-is-not-function
**Symptom:** "Protein kinase" label on catalytically dead pseudokinase. **Fix:** A domain implies a capability, not the substrate or activity. Demote specificity as identity/coverage fall.

### GO IEA circularity in enrichment
**Symptom:** Enrichment results "too clean" — testing the pipeline against itself. **Fix:** Match background to foreground pipeline; prefer non-IEA where possible; pin GO release; caveat results. See `skills/gene-family-evolution/SKILL.md` §7.

### Over-specific transfer
**Symptom:** Exact substrate/EC copied onto distant relative in mechanistically-diverse superfamily. **Fix:** Demote to superfamily / partial EC (`1.1.1.-`); record identity and coverage.

### Wrong eggNOG tax_scope
**Symptom:** Vague or missing function. **Fix:** `--tax_scope auto`, or pin known clade.

### Specialized function from generic hit
**Symptom:** Wrong CAZyme substrate / AMR phenotype / BGC product from plain Pfam. **Fix:** Route to dbCAN / AMRFinderPlus / antiSMASH.

## Quantitative Thresholds

| Threshold | Rationale |
|-----------|-----------|
| Bidirectional coverage ≥ 50–70% query + subject | One-domain coverage = domain-level claim only |
| ~40% identity over full length (well-behaved families) | Soft floor; demote in diverse superfamilies |
| Named fraction > 90% on novel isolate | Over-annotation smell test; expect 20–50% hypothetical |
| eggNOG `--tax_scope auto` | Per-seed informative ceiling |
| KofamScan `*` (adaptive threshold) | Per-KO stricter than global e-value |
| Think bits-per-residue, not raw e-value | e-value scales with DB size; bits/residue = density |

## Verification

- [ ] eggNOG-mapper v2.1.15 pinned; DB version recorded.
- [ ] InterProScan version + member DB release recorded; lookup service or `-dp` noted.
- [ ] Both eggNOG and InterProScan run (orthogonal evidence); results reconciled with provenance.
- [ ] Product names honest: "hypothetical protein" not forced guesses.
- [ ] GO evidence codes tracked; IEA acknowledged.
- [ ] KofamScan adaptive thresholds used for KEGG work.
- [ ] CAZymes / BGCs / AMR routed to specialized tools, not inferred from generic Pfam.

## Related Skills

- `./braker3-pipeline/SKILL.md` — Produces the protein FASTA to annotate
- `./funannotate-pipeline/SKILL.md` — Fungal annotation producing input proteins
- `./annotation-qc/SKILL.md` — Sanity-check annotation coverage before functional work
- `skills/gene-family-evolution/SKILL.md` — §7 uses GO/KEGG from this skill for enrichment of expanded/contracted families
- `skills/ortholog-inference/SKILL.md` — OrthoFinder HOGs contextualize which orthologs are being annotated

## References

- Cantalapiedra CP, et al. 2021. eggNOG-mapper v2. *Mol Biol Evol* 38:5825.
- Huerta-Cepas J, et al. 2019. eggNOG 5.0. *Nucleic Acids Res* 47:D309.
- Jones P, et al. 2014. InterProScan 5. *Bioinformatics* 30:1236.
- Blum M, et al. 2025. InterPro 2025. *Nucleic Acids Res* 53:D444.
- Aramaki T, et al. 2020. KofamKOALA. *Bioinformatics* 36:2251.
- Schnoes AM, et al. 2009. Misannotation in enzyme superfamilies. *PLoS Comput Biol* 5:e1000605.
- Gilks WR, et al. 2002. Percolation of annotation errors. *Bioinformatics* 18:1641.
- Radivojac P, et al. 2013. CAFA. *Nat Methods* 10:221.
- Haynes WA, et al. 2018. Gene annotation bias. *Sci Rep* 8:1362.
