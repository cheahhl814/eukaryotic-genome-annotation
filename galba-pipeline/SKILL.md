---
version: 2
updated: "2026-06-29"
triggers:
  - "galba"
  - "GALBA pipeline"
  - "protein evidence gene prediction"
name: galba-pipeline
description: Run GALBA for protein-only eukaryotic gene prediction using close-relative proteomes. Preferred over BRAKER2 when high-quality close-relative proteins are available.
---

# GALBA Pipeline

**"Predict genes with only protein evidence from close relatives"** -> Use miniprot-aligned close proteomes to train AUGUSTUS, producing more accurate models than BRAKER2 on large vertebrate genomes with good proteome coverage.

## Prerequisites

- **Environment**: Active Pixi environment. Refer to `skills/pixi-skill/SKILL.md`.
- **Upstream evidence**:
  - Soft-masked assembly (`.fasta`). See `./repeat-annotation/SKILL.md`.
  - High-quality proteomes from **close relatives** (`.faa`), not a broad OrthoDB clade.
- **Assembly QC**: BUSCO (genome mode), N50, contamination check passed.

## Installation

```bash
pixi add galba augustus genemark miniprot
```

GALBA bundles miniprot for protein alignment and trains AUGUSTUS on the aligned hints.

## Procedure

### 1. Run GALBA

```bash
galba.pl \
  --genome=masked.fasta \
  --prot_seq=close_relatives.faa \
  --species=my_species \
  --threads=16 \
  --workingdir=galba_out
```

### 2. Key Agentic Flags

| Flag           | Meaning                                             |
| -------------- | --------------------------------------------------- |
| `--genome`     | Soft-masked assembly (mandatory)                    |
| `--prot_seq`   | Close-relative proteomes (FASTA), not broad OrthoDB |
| `--species`    | Unique species name for AUGUSTUS parameter set      |
| `--threads`    | Parallel threads                                    |
| `--workingdir` | Output directory                                    |

### 3. Outputs

| File                       | Content                    |
| -------------------------- | -------------------------- |
| `galba.gtf` / `galba.gff3` | Final gene models          |
| `galba.codingseq`          | Coding sequences           |
| `galba.aa`                 | Predicted proteins (FASTA) |

## When to Prefer GALBA over BRAKER3

| Scenario                  | Use GALBA                     | Use BRAKER3                        |
| ------------------------- | ----------------------------- | ---------------------------------- |
| Evidence                  | Protein only, close relatives | RNA-seq + protein, or RNA-seq only |
| Close proteomes available | Yes (e.g., sister species)    | No / distant only                  |
| Genome size               | Large vertebrate              | Any                                |

**Rule of thumb:** If RNA-seq exists, always prefer BRAKER3 — intron evidence substantially improves splice-site accuracy. GALBA is the best protein-only option when close relatives are available.

## GALBA-Specific Pitfalls

### Using a broad OrthoDB clade instead of close relatives

**Symptom:** Noisier hints, lower accuracy than BRAKER2. **Fix:** GALBA wants a few **close-relative proteomes** (e.g., 3–5 sister species), not a broad clade partition.

### Running on unmasked assembly

**Symptom:** Inflated gene count, spurious models in TEs. **Fix:** Soft-mask first; GALBA does not have a `--softmasking` flag but respects lowercase if the underlying AUGUSTUS is configured correctly. Verify masking upstream.

### Expecting isoforms/UTRs

**Symptom:** mRNA:gene ratio == 1.00. **Fix:** GALBA emits one CDS-only model per locus. Add PASA/Iso-Seq update if isoforms needed.

### Low accuracy on distant relatives

**Symptom:** Systematically mis-placed exons. **Fix:** If close relatives are too distant (>100 Myr), switch to BRAKER2 (broad OrthoDB) or add RNA-seq and run BRAKER3.

## Verification

- [ ] Genome is soft-masked before prediction.
- [ ] Assembly QC passed before prediction.
- [ ] Input proteomes are from close relatives (not broad OrthoDB).
- [ ] GALBA completed without errors.
- [ ] Output sanity-checked (gene count, mono-exonic fraction, mRNA:gene ratio).
- [ ] Protein-BUSCO vs genome-BUSCO comparison shows no major gene drop.
- [ ] If isoforms/UTRs needed, PASA/Iso-Seq update step is planned.

## Related Skills

- `./gene-prediction-decision-tree/SKILL.md` — Router: choose the right tool
- `./repeat-annotation/SKILL.md` — PREREQUISITE: soft-mask assembly
- `./braker3-pipeline/SKILL.md` — Default when RNA-seq + protein exist
- `./annotation-qc/SKILL.md` — Assess annotation quality post-prediction

## References

- Brůna T, et al. 2023. GALBA: genome annotation with miniprot and AUGUSTUS. *BMC Bioinformatics* 24:327.
- Gabriel L, et al. 2024. BRAKER3: fully automated genome annotation using RNA-seq and protein evidence. *Genome Res* 34:769-777.
