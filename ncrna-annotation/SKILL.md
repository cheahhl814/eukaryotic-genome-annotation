---
version: 2
updated: "2026-06-29"
triggers:
  - "ncRNA annotation"
  - "tRNA rRNA annotation"
  - "infernial rfam"
  - "censuscan rfam"
name: ncrna-annotation
description: Identify non-coding RNAs (tRNA, rRNA, snoRNA, snRNA, miRNA, lncRNA) using Infernal/Rfam covariance models, tRNAscan-SE 2.0, barrnap, and expression-based methods. Covers structure-aware homology search and why BLAST fails for ncRNA.
---

# Non-Coding RNA Annotation

**"Find non-coding RNAs in my genome"** -> Scan an assembly for structured ncRNA families using covariance models (sequence + secondary structure jointly), with specialist detectors for tRNA, rRNA, and expression-based boundaries for miRNA/lncRNA.

## The Single Most Important Modern Insight

ncRNA homology is **structure-conserved, not sequence-conserved**. A G-C pair can become A-U across evolution while the *structure* is preserved — both positions mutate together (compensatory substitution). To a sequence-only tool these look like two mismatches; to a covariance model (the Rfam/Infernal engine) the correlated change is the strongest possible evidence of homology.

**Three consequences:**
1. **BLAST is categorically wrong for structured ncRNA.** Two RNase P RNAs may share <60% identity (BLAST sees noise) yet have unmistakable shared structure. Use Infernal covariance models, not BLAST.
2. **Homology-based annotation is a recall floor, never a count.** Fast-evolving and lineage-specific ncRNAs are invisible by design. Report "at least N conserved-family loci," never "the genome has N ncRNAs."
3. **Whole classes need expression evidence, not genomic search.** miRNAs need small-RNA-seq; lncRNAs are transcript catalogs, not structurally annotatable.

## Tool Taxonomy

| RNA class | Tool | Method |
| --- | --- | --- |
| **tRNA** | tRNAscan-SE 2.0 | Isotype-specific CMs + pseudogene logic |
| **rRNA (fast)** | barrnap | nhmmer HMM profiles; kingdom flag |
| **rRNA (structure-aware)** | Infernal + Rfam SSU/LSU | CM search; better boundaries / unusual taxa |
| **tmRNA** | ARAGORN | Heuristic cloverleaf + tmRNA models |
| **miRNA** | miRDeep2 (+ small-RNA-seq) | Dicer-processing model on read pileups |
| **snoRNA** | snoscan / snoReport | Guide-target complementarity / SVM |
| **Everything else structured** | Infernal `cmscan` vs Rfam | Covariance models; the general engine |
| **lncRNA** | StringTie + CPC2/CPAT/FEELnc | Transcript assembly + coding-potential |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
| --- | --- | --- |
| Prokaryote, fast complete annotation | barrnap + tRNAscan-SE + Infernal/Rfam | What Bakta/Prokka wrap |
| tRNA specifically | **tRNAscan-SE 2.0** | Isotype, pseudogene, intron logic; not Rfam's generic tRNA model |
| rRNA, speed matters | **barrnap** | Seconds per genome |
| Broad ncRNA sweep | **Infernal `cmscan`** vs full Rfam (GA + clan competition) | Structure-aware, family-typed |
| miRNA | **miRDeep2 + small-RNA-seq** | Genomic hairpin prediction is unreliable |
| lncRNA | **Transcript assembly + coding-potential** | Not structurally conserved; no CM |
| Claim a conserved structure | **R-scape** covariation test | Thermodynamic fold ≠ selected structure |

## Prerequisites

- **Environment**: Active Pixi environment. Refer to `skills/pixi-env-mgmt/SKILL.md`.
- **Input**: Assembled genome (`.fa`) or transcriptome for lncRNA.
- **Databases**: Rfam covariance model library (`.cm`), pre-pressed with `cmpress`.

## Installation

```bash
pixi add infernal trnascan-se barrnap aragorn
# For expression-based miRNA/lncRNA:
pixi add mirdeep2 stringtie gffread
# Python parsing:
pixi add pandas
```

Download and prepare Rfam (do once per release):
```bash
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm
# Also download clan information for deoverlap:
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
```

## Procedure

### 1. Infernal / cmscan — The General ncRNA Engine

```bash
# Full canonical command with all correctness flags
cmscan -Z $(echo "$(grep -v '>' genome.fa | tr -d '\n' | wc -c) * 2 / 1000000" | bc) \
       --cut_ga --rfam --nohmmonly \
       --tblout out.tblout --fmt 2 --clanin Rfam.clanin \
       Rfam.cm genome.fa > out.cmscan

# Drop within-clan overlaps (mandatory — one locus hits several related models)
grep -v " = " out.tblout > out.deoverlapped.tblout
```

**Key flags:**
| Flag | Meaning |
|------|---------|
| `--cut_ga` | Per-family gathering threshold (curator-set bit-score floor). **Never override** — a fixed E-value treats 70-nt tRNA and 2,900-nt LSU identically. |
| `--nohmmonly` | Force full CM (structure-aware) scoring so scores are GA-comparable. |
| `-Z` | Database size in Mb for run-comparable E-values: `residues × 2 / 1e6` (both strands). |
| `--fmt 2 --clanin` | Enable clan competition + `grep -v " = "` deoverlap. **Mandatory** — raw table double-counts within clans (one 16S locus becomes "several rRNA genes"). |

### 2. tRNAscan-SE 2.0 — The tRNA Authority

```bash
# Bacterial mode (use -E for eukaryote, -A for archaea, -G for mixed/metagenome)
tRNAscan-SE -B -o trna.out -f trna.ss -m trna.stats --gff trna.gff3 genome.fa
```

**Domain modes:**
| Flag | Domain |
|------|--------|
| `-E` | Eukaryotic (default) |
| `-B` | Bacterial |
| `-A` | Archaeal |
| `-G` | General / metagenome |
| `-O` | Organelar (disables pseudogene checking) |

> **Critical:** Domain choice is **not cosmetic** — wrong mode mis-scores isotypes; there is no auto-detect. Report the **high-confidence set**, not raw hits (raw counts include pseudogenes/SINEs and can be 2–10× inflated in eukaryotes).

### 3. barrnap — Fast rRNA Detection

```bash
barrnap --kingdom euk genome.fa > rrna.gff3   # euk | bac | arc | mito
```

> **rDNA copy number from an assembly is essentially always wrong** — near-identical arrays collapse in short-read assemblies. The annotated count is a floor (off by orders of magnitude in eukaryotes); use long reads or read depth for true copy number.

### 4. Combining Results with Python

Merge Infernal and tRNAscan-SE into one annotation, **preferring the specialist for tRNA** (drop Infernal tRNA rows):

```python
import pandas as pd

def parse_infernal_tbl(tbl_file):
    rows = []
    with open(tbl_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 18:
                continue
            rows.append({
                'rfam_name': p[1], 'seqid': p[3], 'strand': p[11],
                'score': float(p[16]), 'evalue': float(p[17])
            })
    df = pd.DataFrame(rows)
    # tRNAscan-SE owns tRNA — drop Rfam generic tRNA hits
    return df[~df['rfam_name'].str.contains('tRNA', case=False, na=False)]
```

## Pitfalls

### BLAST for structured ncRNA
**Symptom:** Misses structure-conserved homologs; ranks pseudogenes above real homologs. **Fix:** Always use Infernal/Rfam covariance models.

### Missing `--cut_ga` or clan competition
**Symptom:** Inflated/redundant family calls; one locus counted several times. **Fix:** Full canonical command + `grep -v " = "` deoverlap.

### Wrong tRNAscan-SE domain mode
**Symptom:** Mis-scored/miscalled isotypes. **Fix:** `-B` / `-A` / `-G` per source organism; never use default `-E` on prokaryotes.

### Counting raw tRNA/rRNA hits
**Symptom:** 150+ bacterial "tRNAs"; "5 copies of 18S" in eukaryotes. **Fix:** Report tRNAscan-SE **high-confidence set**; treat rRNA counts as floors; flag NUMT-type organellar tRNAs.

### miRNA from genome alone
**Symptom:** Thousands of false-positive "miRNAs." **Fix:** Demand small-RNA-seq Dicer signature (miRDeep2). Label homology-only as "candidates."

### lncRNA from CM search
**Symptom:** No results or garbage. **Fix:** lncRNAs have no conserved structure — use transcript assembly + coding-potential (CPC2/CPAT/FEELnc).

## Quantitative Thresholds

| Threshold | Rationale |
| --- | --- |
| `--cut_ga` per-family GA threshold | Family-specific noise floor; never override genome-wide |
| `-Z` = residues × 2 / 1e6 | Run-comparable E-values (both strands) |
| Bacterial tRNA ~28–90 (E. coli ~89) | 150+ implies fragments/pseudogenes |
| Eukaryotic tRNA ~170–570 (high-confidence) | Raw hits far higher (SINEs/pseudogenes) |
| Bacterial rRNA operons 1–15 | Annotated count is a floor (rDNA collapse) |
| miRNA needs small-RNA-seq Dicer signature | Foldability is not a miRNA; >1/3 miRBase is artifact |
| R-scape before claiming conserved structure | Thermodynamic fold ≠ selected structure |

## Verification

- [ ] `--cut_ga` used in Infernal (not fixed E-value).
- [ ] `--fmt 2 --clanin` + `grep -v " = "` deoverlap applied.
- [ ] tRNAscan-SE domain mode matches organism (`-B`/`-E`/`-A`/`-G`).
- [ ] High-confidence tRNA set reported (not raw hits).
- [ ] rRNA copy numbers treated as floors (assembly-collapse aware).
- [ ] miRNA claims backed by small-RNA-seq Dicer evidence.
- [ ] lncRNAs from transcript assembly, not CM search.
- [ ] Rfam version recorded (results vary by release).

## Related Skills

- `./repeat-annotation/SKILL.md` — TEs and ncRNAs are distinct; don't confuse them
- `./gene-prediction-decision-tree/SKILL.md` — Protein-coding predictors miss ncRNAs
- `./braker3-pipeline/SKILL.md` — BRAKER does not find ncRNAs
- `./annotation-qc/SKILL.md` — tRNA/rRNA count sanity in QC panel
- `skills/read-qc-trimming/SKILL.md` — Small-RNA-seq quality for miRNA work

## References

- Nawrocki EP, Eddy SR. 2013. Infernal 1.1: 100-fold faster RNA homology searches. *Bioinformatics* 29:2933–2935.
- Ontiveros-Palacios N, et al. 2025. Rfam 15: RNA families database. *Nucleic Acids Res* 53:D258-D267.
- Chan PP, et al. 2021. tRNAscan-SE 2.0. *Nucleic Acids Res* 49:9077–9096.
- Laslett D, Canback B. 2004. ARAGORN, a program to detect tRNA and tmRNA genes. *Nucleic Acids Res* 32:11–16.
- Friedländer MR, et al. 2012. miRDeep2 accurately identifies known and novel microRNA genes. *Nucleic Acids Res* 40:37–52.
- Fromm B, et al. 2015. A uniform system for the annotation of vertebrate microRNA genes (MirGeneDB). *Annu Rev Genet* 49:213–242.
- Rivas E, Clements J, Eddy SR. 2017. A statistical test for conserved RNA structure shows lack of evidence for structure in lncRNAs (R-scape). *Nat Methods* 14:45–48.
- Graur D, et al. 2013. On the immortality of television sets: "function" in the human genome according to ENCODE. *Genome Biol Evol* 5:578–590.
