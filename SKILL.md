---
name: "eukaryotic-genome-annotation"
description: "End-to-end eukaryotic genome annotation (nuclear genomes of fungi, plants, and animals; organelle genomes are out of scope — see mitotree-agent): soft-mask repeats, predict protein-coding genes (BRAKER3/GALBA/Funannotate), annotate ncRNAs, assign functional labels, and QC the result. Routes through a decision tree to select the optimal tool per evidence type and clade; ships a Nextflow DSL2 runner (runners/nextflow-runner/) for reproducible unattended execution."
version: 2
updated: "2026-06-29"
triggers:
  - "annotate genome"
  - "genome annotation"
  - "predict genes"
  - "BRAKER3"
  - "Funannotate"
  - "GALBA"
  - "repeat annotation"
  - "functional annotation"
  - "ncRNA annotation"
  - "gene prediction"
  - "eukaryotic genome annotation"
  - "annotate my genome"
  - "soft-mask repeats"
  - "GO KEGG Pfam annotation"
  - "annotation QC"
requires:
  - "Pixi environment (env managed via this skill's pixi.toml; use skills/pixi-skill for env management procedures)"
  - "Soft-masked genome assembly (produced by Phase 1 if not already done)"
  - "RNA-seq reads and/or close-relative proteomes (optional but strongly recommended)"
  - "biodb-fetch skill (https://github.com/cheahhl814/biodb-fetch — OPTIONAL soft dep: fetches close-relative RefSeq proteomes for protein-evidence gene prediction; degrade gracefully to user-supplied .faa when not deployed/offline)"
---

# Genome Annotation Pipeline

## Agent Compatibility

This skill bundle is agent-agnostic. All procedures are written as Markdown instructions that any LLM coding agent can follow. The only adaptation needed is mapping Pixi/bash tool calls to your agent's native equivalents.

| Tool in this skill               | Pi     | Claude Code | Codex  | Generic Bash |
| -------------------------------- | ------ | ----------- | ------ | ------------ |
| `pixi`                           | ✅ bash | ✅ bash      | ✅ bash | ✅ bash       |
| `braker.pl`                      | ✅ bash | ✅ bash      | ✅ bash | ✅ bash       |
| `galba.pl`                       | ✅ bash | ✅ bash      | ✅ bash | ✅ bash       |
| `funannotate`                    | ✅ bash | ✅ bash      | ✅ bash | ✅ bash       |
| `RepeatModeler/Masker`           | ✅ bash | ✅ bash      | ✅ bash | ✅ bash       |
| `cmscan` / `tRNAscan-SE`         | ✅ bash | ✅ bash      | ✅ bash | ✅ bash       |
| `emapper.py` / `interproscan.sh` | ✅ bash | ✅ bash      | ✅ bash | ✅ bash       |
| `busco` / `compleasm` / `OMArk`  | ✅ bash | ✅ bash      | ✅ bash | ✅ bash       |

All sub-skills are linked via relative paths (e.g., `./repeat-annotation/SKILL.md`). Keep the directory structure intact when importing.

A 6-phase pipeline that takes you from "I have a genome assembly" to a fully annotated genome with protein-coding genes, ncRNAs, functional labels, and QC metrics.

- **Phase 0:** Assembly QC (prerequisite gate)
- **Phase 1:** Repeat Annotation & Soft-Masking
- **Phase 2:** Gene Prediction (routed by evidence type and clade)
- **Phase 3:** ncRNA Annotation
- **Phase 4:** Functional Annotation
- **Phase 5:** Annotation QC

**Phase 0** is a gate — do not proceed if the assembly fails BUSCO/N50/contamination checks. **Phase 1** is mandatory; unmasked genomes produce 2× inflated gene counts. **Phase 2** routes through the decision tree based on your evidence. **Phase 3** runs in parallel with or after Phase 2. **Phase 4** adds meaning to predicted proteins. **Phase 5** decides whether the annotation is publishable or needs rework.

## When to Use

- "Annotate my eukaryotic genome" (the primary use case)
- "My gene count is 2× what I expected" (likely missing Phase 1 or wrong predictor)
- "Which gene predictor should I use?" (→ Phase 2 decision tree)
- "Add GO/KEGG/Pfam to my predicted proteins" (→ Phase 4)
- "Is my annotation any good?" (→ Phase 5)
- Any eukaryotic genome annotation task that needs reproducible, QC'd output

## Prerequisites

| Requirement                       | Install                             | Phase      |
| --------------------------------- | ----------------------------------- | ---------- |
| Pixi environment                  | See `skills/pixi-skill/SKILL.md` | All        |
| Genome assembly (`.fa`)           | User-supplied                       | Phase 0    |
| RNA-seq reads (`.fq.gz`)          | Optional but recommended            | Phase 2    |
| Close-relative proteomes (`.faa`) | Optional; via [biodb-fetch](https://github.com/cheahhl814/biodb-fetch) (NCBI Datasets → RefSeq proteomes) or user-supplied | Phase 2    |
| OrthoDB clade partition           | Download from OrthoDB (bulk partition FASTA — biodb-fetch does NOT cover this; its orthologs tooling is per-gene lookups) | Phase 2    |
| Rfam CM library                   | Download from Rfam (biodb-fetch does not cover Rfam; grab `Rfam.cm` + `Rfam.clanin` from the EBI FTP) | Phase 3    |
| eggNOG DB (~44 GB)                | `download_eggnog_data.py` (tool-canonical; biodb-fetch deliberately out of scope) | Phase 4    |
| InterProScan DB                   | Download separately (tool-canonical) | Phase 4    |
| BUSCO/compleasm lineage dataset   | `odb10` / `odb12` (busco's own `--download`) | Phase 0, 5 |

> **Fallback paths:** Phases 2–5 each have evidence-variant paths. If you lack RNA-seq, use protein-only or DL modes. If you lack both, Helixer/Tiberius can run ab initio on well-represented clades (see decision tree).
>
> **Database acquisition routing (who downloads what):** most annotation DBs are **tool-canonical
> downloads** — eggNOG via `download_eggnog_data.py`, BUSCO lineages via busco's `--download`,
> InterProScan via its setup script, Rfam via the EBI FTP (the [biodb-fetch](https://github.com/cheahhl814/biodb-fetch)
> skill's 8 retrieval domains do not cover these bulk tool DBs, and its per-gene ortholog REST
> tooling is not a substitute for the bulk OrthoDB clade partition BRAKER3 needs). The one fit:
> **close-relative proteomes** — fetch RefSeq reference proteomes of sister species through
> biodb-fetch's NCBI Datasets wrapper (md5-verified, store-cached, `manifest.tsv` contract v2);
> see `galba-pipeline/SKILL.md` §Prerequisites for the command. Optional route — degrade to
> user-supplied `.faa` when biodb-fetch is unavailable.

## Procedure

### PHASE 0 — ASSEMBLY QC (Gate)

**Do not skip.** Annotating a bad assembly wastes compute and produces confidently-wrong gene models.

1. Run BUSCO/compleasm in **genome mode** on the assembly:
   
   ```bash
   busco -i assembly.fa -m genome -l vertebrata_odb12 -o busco_asm -c 16
   # or faster:
   compleasm run -a assembly.fa -l vertebrata -o compleasm_asm -t 16
   ```
2. Check N50, contamination (FCS-GX / BlobTools).
3. **Gate criteria:**
   - BUSCO completeness ≥ 90% (adjust per clade)
   - Contamination < 5%
   - N50 appropriate for clade
4. If the assembly fails, **stop and fix the assembly** before proceeding.

### PHASE 1 — REPEAT ANNOTATION & SOFT-MASKING

**Mandatory prerequisite for all gene prediction.** Unmasked TEs inflate gene counts 2× and produce thousands of spurious single-exon models.

→ See [repeat-annotation](./repeat-annotation/SKILL.md) for the full procedure.

**Quick path (non-model default):**

```bash
# EarlGrey wraps RepeatModeler2 + curation + masking
earlgrey -g assembly.fa -o earlgrey_out -t 16 -s my_species
```

**Key output:** `assembly.fa.masked` (soft-masked genome, lowercase = repeats)

> ⚠️ **Critical:** Always use soft-masking (`-xsmall`). Hard-masking (`N`) silently truncates genes overlapping repeats. Decontaminate the repeat library against a protein DB before masking to avoid the "domesticated-gene massacre."

### PHASE 2 — GENE PREDICTION

Route through the decision tree based on your evidence and clade.

→ See [gene-prediction-decision-tree](./gene-prediction-decision-tree/SKILL.md) for the full router.

| Evidence / Scenario                    | Tool                            | Sub-skill                                               |
| -------------------------------------- | ------------------------------- | ------------------------------------------------------- |
| RNA-seq + protein DB                   | **BRAKER3** (full evidence)     | [braker3-pipeline](./braker3-pipeline/SKILL.md) §2      |
| RNA-seq only                           | **BRAKER3** (RNA-seq-only mode) | [braker3-pipeline](./braker3-pipeline/SKILL.md) §2c     |
| Protein only, close relatives          | **GALBA**                       | [galba-pipeline](./galba-pipeline/SKILL.md)             |
| Protein only, broad/distant            | **BRAKER3** (protein-only mode) | [braker3-pipeline](./braker3-pipeline/SKILL.md) §2b     |
| Fungus (any evidence)                  | **Funannotate**                 | [funannotate-pipeline](./funannotate-pipeline/SKILL.md) |
| No evidence, represented clade         | Helixer / Tiberius              | Decision-tree only                                      |
| Reference annotation of close relative | GeMoMa                          | Decision-tree only                                      |

**Key output:** `braker.gff3` / `galba.gff3` / `funannotate.gff3` + `proteins.faa`

> ⚠️ **Pitfall:** BRAKER3, GALBA, and Funannotate all emit one CDS-only model per locus (mRNA:gene = 1.00). If isoforms/UTRs are needed, add a PASA/Iso-Seq update step.

### PHASE 3 — ncRNA ANNOTATION

Runs in parallel with or after Phase 2. Protein-coding predictors do not find ncRNAs.

→ See [ncrna-annotation](./ncrna-annotation/SKILL.md) for the full procedure.

**Minimum path:**

```bash
# Infernal (general structured ncRNAs)
cmscan --cut_ga --rfam --nohmmonly --fmt 2 --clanin Rfam.clanin Rfam.cm genome.fa > ncRNA.cmscan

# tRNAscan-SE (the tRNA authority)
tRNAscan-SE -E -o trna.out --gff trna.gff3 genome.fa

# barrnap (fast rRNA)
barrnap --kingdom euk genome.fa > rrna.gff3
```

**Key output:** Combined GFF3 of tRNA, rRNA, snoRNA, snRNA, miRNA annotations.

> ⚠️ **Critical:** Use `--cut_ga` (never fixed E-value). Apply clan competition deoverlap (`grep -v " = "`). Report tRNAscan-SE high-confidence set, not raw hits. miRNAs require small-RNA-seq evidence.

### PHASE 4 — FUNCTIONAL ANNOTATION

Assign GO terms, Pfam/InterPro domains, KEGG orthologs, EC numbers, and product names to predicted proteins.

→ See [functional-annotation](./functional-annotation/SKILL.md) for the full procedure.

**Default pair:** eggNOG-mapper (orthology) + InterProScan (domain signatures) — orthogonal evidence that dampens single-source error.

```bash
# eggNOG-mapper
emapper.py -i proteins.faa --itype proteins -m diamond --tax_scope auto --cpu 16 -o eggnog_out

# InterProScan
interproscan.sh -i proteins.faa -f TSV,GFF3 -goterms -pa -cpu 16
```

**Key output:** Annotated protein FASTA + GAF (GO evidence codes) + eggNOG annotations + InterPro entries.

> ⚠️ **Pitfall:** "Domain presence ≠ function known." A Pfam "protein kinase" domain on a catalytically dead pseudokinase. Record provenance. Route CAZymes/BGCs/AMR to specialized tools (dbCAN/antiSMASH/AMRFinderPlus).

### PHASE 5 — ANNOTATION QC

The single most important check most people skip: **run BUSCO/compleasm in both genome and protein mode on the same assembly.**

→ See [annotation-qc](./annotation-qc/SKILL.md) for the full procedure.

**The diagnostic comparison:**

| Scenario               | Assembly BUSCO | Proteome BUSCO | Meaning                          |
| ---------------------- | -------------- | -------------- | -------------------------------- |
| Good annotation        | 98%            | 95%+           | ✅ Genes present and found        |
| Predictor missed genes | 98%            | 85%            | ⚠️ Fix masking/evidence/training |
| Bad assembly           | 80%            | 78%            | ❌ Stop; fix the assembly first   |

```bash
busco -i proteins.faa -m proteins -l vertebrata_odb12 -o busco_prot -c 16
busco -i assembly.fa -m genome -l vertebrata_odb12 -o busco_genome -c 16
```

Add OMArk for proteome consistency (catches over-prediction BUSCO misses) and the gene-set sanity panel (mono-exonic fraction, protein-length distribution, mRNA:gene ratio).

## Pitfalls

- **Skipping Phase 0:** Annotating a contaminated or fragmented assembly produces confidently-wrong gene models. Always QC first.
- **Skipping Phase 1:** Unmasked TEs cause 2× inflated gene counts and thousands of spurious single-exon genes. Soft-masking is mandatory.
- **Hard-masking:** Using `N` masking silently truncates genes overlapping repeats. Always use `-xsmall`.
- **Wrong OrthoDB partition:** Systematically mis-placed exon boundaries with clean-looking GFF3. Pick the smallest clade partition containing your species.
- **Expecting isoforms from a one-model predictor:** BRAKER3/GALBA emit one CDS-only model per locus. Add PASA/Iso-Seq for UTRs and isoforms.
- **High BUSCO-Duplicated celebrated as success:** D > 5–8% without known WGD = uncollapsed haplotigs. Run `purge_dups` before annotating.
- **BLAST for ncRNA:** BLAST misses structure-conserved homologs. Always use Infernal/Rfam covariance models.
- **Functional over-annotation:** 95% "named" proteins on a novel genome is a red flag. Expect 20–50% hypothetical. Record provenance and evidence codes.
- **GO IEA circularity in enrichment:** Testing the annotation pipeline against itself. Match background to foreground pipeline; prefer non-IEA.
- **Counting genes as a quality metric:** Gene count varies 20–40% between pipelines on the same assembly. The signal is in mono-exonic fraction, protein-length distribution, and mRNA:gene ratio.

## Verification

1. Assembly passed Phase 0 QC (BUSCO ≥ 90%, contamination < 5%).
2. Genome is soft-masked (lowercase letters in FASTA, no `N` masking).
3. Repeat library decontaminated against protein DB (no domesticated-gene massacre).
4. Correct gene predictor selected by evidence type and clade (Phase 2 decision tree).
5. ncRNA annotation uses `--cut_ga` with clan competition deoverlap (Phase 3).
6. Both genome-mode and proteome-mode BUSCO/compleasm run (Phase 5 diagnostic).
7. OMArk consistency check passed (eukaryotes).
8. Gene-set sanity panel computed (mono-exonic %, protein length, mRNA:gene).
9. Functional annotations carry provenance (eggNOG version, InterProScan version, evidence codes).
10. CAZymes/BGCs/AMR routed to specialized tools, not inferred from generic Pfam.

## Workflow vs Manual Execution

The agentic skills (this document and sub-skills) handle **judgment, routing, and interpretation**. For **reproducible, unattended execution** of Phases 1–5, use the Nextflow workflow:

→ See [runners/nextflow-runner/SKILL.md](./runners/nextflow-runner/SKILL.md) for when to use the workflow vs manual execution.

```bash
# After Phase 0 assembly QC gate passes:
nextflow run runners/nextflow-runner/main.nf \
  -profile singularity,slurm \
  --genome assembly.fa \
  --mode braker3 \
  --species 'My_species' \
  --orthodb OrthoDB/Vertebrata.fa \
  --rnaseq 'reads/*_{1,2}.fq.gz' \
  --outdir results/
```

**What the workflow handles:** deterministic execution, checkpointing, resume, containerization, resource management.
**What the agent handles:** Phase 0 QC gate, tool selection, OrthoDB partition choice, repeat library curation, QC interpretation, failure diagnosis.

## Related Skills

- [repeat-annotation](./repeat-annotation/SKILL.md) — Phase 1: soft-mask repeats
- [gene-prediction-decision-tree](./gene-prediction-decision-tree/SKILL.md) — Phase 2: router
- [braker3-pipeline](./braker3-pipeline/SKILL.md) — Phase 2: default gene predictor
- [galba-pipeline](./galba-pipeline/SKILL.md) — Phase 2: protein-only close relatives
- [funannotate-pipeline](./funannotate-pipeline/SKILL.md) — Phase 2: fungal standard
- [ncrna-annotation](./ncrna-annotation/SKILL.md) — Phase 3: non-coding RNAs
- [functional-annotation](./functional-annotation/SKILL.md) — Phase 4: GO/KEGG/Pfam/EC
- [annotation-qc](./annotation-qc/SKILL.md) — Phase 5: quality control
- [runners/nextflow-runner/SKILL.md](./runners/nextflow-runner/SKILL.md) — Nextflow workflow for reproducible execution
- `skills/pixi-skill/SKILL.md` — Environment management
- `skills/rna-seq-skills/read-mapping-alignment/SKILL.md` — Splice-aware RNA-seq alignment (upstream)
- biodb-fetch (https://github.com/cheahhl814/biodb-fetch) — OPTIONAL: close-relative RefSeq proteome retrieval (protein-evidence Phase 2); not a hard dependency
- `skills/nf-skill/SKILL.md` — Nextflow DSL2 conventions and patterns
- ortholog-inference (planned, not yet in @skills/) — OrthoFinder HOGs (downstream)
- gene-family-evolution (planned, not yet in @skills/) — GO enrichment on expanded/contracted families (downstream)