---
version: 2
updated: "2026-06-29"
triggers:
  - "funannotate"
  - "fungal annotation"
  - "fungal genome annotation"
name: funannotate-pipeline
description: Run Funannotate for fungal genome annotation. Covers masking, training, prediction, and PASA-based isoform/UTR updates — the fungal de facto standard.
---

# Funannotate Pipeline

**"Annotate my fungal genome"** -> Use the fungal de facto standard pipeline: mask -> train -> predict (EVM) -> update (PASA) for gene models with UTRs and isoforms.

## Prerequisites

- **Environment**: Active Pixi environment. Refer to `skills/pixi-env-mgmt/SKILL.md`.
- **Upstream evidence**:
  - Assembled genome (`.fa`). Funannotate handles masking internally but soft-masking is preferred upstream.
  - RNA-seq reads (`.fq.gz`) or transcript evidence (`.fa`) — optional but strongly recommended.
  - Protein evidence from related fungi (`.fa`) — optional.
- **Assembly QC**: Fungal genomes are smaller; check for contamination and completeness.

## Installation

### Via Pixi

```bash
pixi add funannotate hisat2 samtools augustus genemark busco
# Funannotate bundles EVM, PASA, and fungal-specific training logic
```

### Via Docker (Alternative)

```bash
# Funannotate provides official containers; use if dependency resolution fails
docker pull nextgenusfs/funannotate:latest
```

## Procedure

### 1. Mask Repeats

Funannotate can mask internally, but preferred to use a curated soft-masked assembly from `./repeat-annotation/SKILL.md`.

```bash
# Or let Funannotate mask (if no prior soft-mask exists)
funannotate mask -i assembly.fa -o masked.fa --cpus 16
```

### 2. Train (RNA-seq Evidence)

```bash
funannotate train \
  -i masked.fa \
  -o out_dir \
  -l R1.fq.gz -r R2.fq.gz \
  --species "Genus species" \
  --cpus 16
```

### 3. Predict Genes

```bash
funannotate predict \
  -i masked.fa \
  -o out_dir \
  -s "Genus species" \
  --transcript_evidence transcripts.fa \
  --protein_evidence proteins.fa \
  --cpus 16
```

### 4. Update with UTRs and Isoforms (PASA)

```bash
funannotate update \
  -i out_dir \
  --cpus 16
```

### Key Agentic Flags

| Step | Flag | Meaning |
|------|------|---------|
| `mask` | `-i` | Input assembly |
| `train` | `-l / -r` | Left/right RNA-seq reads |
| `predict` | `--transcript_evidence` | Transcript FASTA (training) |
| `predict` | `--protein_evidence` | Protein FASTA (hinting) |
| `update` | PASA-based | Adds UTRs and isoforms |

## Fungal-Specific Pitfalls

### Tiny introns
Fungi have introns as short as ~30–40 bp. Standard vertebrate-trained predictors miss these. Funannotate is **tiny-intron-aware** — do not substitute generic BRAKER3 for fungal annotation.

### Over-masking with plant/animal repeat libraries
**Symptom:** Real fungal genes lost to over-aggressive masking. **Fix:** Use a fungal-specific repeat library or Dfam Fungi partition.

### Missing RNA-seq leading to poor training
**Symptom:** Low BUSCO completeness. **Fix:** Even a small RNA-seq dataset dramatically improves training. If no RNA-seq, use `--protein_evidence` from closely related fungi.

### Skipping the `update` step
**Symptom:** mRNA:gene == 1.00, no UTRs, isoform analyses fail. **Fix:** Always run `funannotate update` if downstream analyses need UTRs or isoforms.

### Expecting gene names from de novo prediction
**Symptom:** Output has generic IDs (`FUN_000001`). **Fix:** Funannotate assigns standard IDs; functional naming requires `funannotate annotate` with InterProScan/eggNOG. See `./functional-annotation/SKILL.md`.

## Verification

- [ ] Assembly is masked (internally or upstream).
- [ ] Fungal-specific repeat library used (not plant/animal).
- [ ] `train` step completed with evidence (RNA-seq or proteins).
- [ ] `predict` step produced GFF3 with expected gene count (~5,000–15,000 for typical fungi).
- [ ] `update` step run if UTRs/isoforms needed.
- [ ] Protein-BUSCO score matches or exceeds genome-BUSCO.
- [ ] If functional naming needed, `funannotate annotate` or equivalent planned.

## Related Skills

- `./gene-prediction-decision-tree/SKILL.md` — Router: choose the right tool
- `./repeat-annotation/SKILL.md` — Soft-mask assembly before prediction
- `./functional-annotation/SKILL.md` — Add GO/KEGG/Pfam to predicted proteins
- `./annotation-qc/SKILL.md` — Assess annotation quality

## References

- Palmer JM, Stajich J. 2020. Funannotate v1.8: eukaryotic genome annotation. Zenodo. doi:10.5281/zenodo.4054262.
- Haas BJ, et al. 2008. Automated eukaryotic gene structure annotation using EVidenceModeler. *Genome Biol* 9:R7.
- Haas BJ, et al. 2003. Improving the Arabidopsis genome annotation using maximal transcript alignment assemblies (PASA). *Nucleic Acids Res* 31:5654-5666.
