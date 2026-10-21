# Genome Annotation Pipeline (Eukaryotic)

[![Version](https://img.shields.io/badge/version-1.0.0-blue)](#-installation)
[![Type](https://img.shields.io/badge/type-agent%20skill-blueviolet)](#-installation)
[![Built with](https://img.shields.io/badge/built%20with-bioinfo--skill--creator-orange)](https://github.com/cheahhl814/skill-forge)

End-to-end eukaryotic genome annotation pipeline (nuclear genomes of fungi, plants, and animals; organelle genomes are out of scope — see mitotree-agent). Soft-mask repeats, predict protein-coding genes (BRAKER3/GALBA/Funannotate), annotate ncRNAs, assign functional labels, and QC the result. Routes through a decision tree to select the optimal tool per evidence type and clade, and ships a Nextflow DSL2 runner (runners/nextflow-runner/) for reproducible unattended execution.

**Repository**: https://github.com/cheahhl814/eukaryotic-genome-annotation

> [!NOTE]
> Current version: **v1.0.0** (updated 2026-10-21).

## Contents

- [Installation](#-installation)
- [Usage](#-usage)
- [Pipeline overview](#-pipeline-overview)
- [Tools](#-tools)
- [Update check](#-update-check)
- [Repository layout](#-repository-layout)
- [Hard guarantees](#-hard-guarantees)
- [Nextflow runner](#-nextflow-runner)
- [Provenance](#-provenance)

## 🚀 Installation

This is an **agent skill**, not a user-facing library. The recommended install path is to let your AI agent import it.

**Option A — give your agent this prompt (recommended):**

```text
Install the eukaryotic-genome-annotation skill from
https://github.com/cheahhl814/eukaryotic-genome-annotation —
clone it into your agent's skills directory (the path your agent watches
for skills) and run `pixi install` from that directory. Then read the
skill's SKILL.md to understand its phases and the Input/Output contract
for each. Confirm when the environment is ready.
```

**Option B — manual install:**

```bash
git clone https://github.com/cheahhl814/eukaryotic-genome-annotation.git
cd eukaryotic-genome-annotation
pixi install
```
> [!TIP]
> `pixi install` resolves every pinned tool from `pixi.toml` (channels: conda-forge, bioconda) into an isolated `.pixi/` environment — no system-wide installs, no version conflicts with other skills.

## 💡 Usage

The skill is designed to be driven by an AI agent: the agent reads the master `SKILL.md`, detects the current phase from the filesystem state of your run directory, and routes to the right sub-skill.

### Natural-language prompts that trigger the skill

```text
annotate genome
genome annotation
predict genes
BRAKER3
Funannotate
```
### Manual phase execution

If you prefer to drive the phases yourself (or want to re-run a single phase):

```bash
# Set your run directory (all phase artifacts land here)
RUN_DIR=/path/to/run-dir

pixi run preflight   # Phase 1 — validate inputs, write preflight.md + params.json
pixi run run         # Phase 2 — execute the workflow (after preflight ≥ GO)
pixi run qc          # Phase 3 — build the final report (after run completes)
pixi run debug       # On failure — interpret stderr via the signature library
pixi run battle-test # Verify the skill's structural integrity
```

> [!IMPORTANT]
> Each phase has an explicit **Inputs/Outputs contract** at the top of its sub-skill `SKILL.md`. If the upstream artifact is missing (e.g. you run `run` before `preflight` passed), the sub-skill refuses to proceed and tells you which phase to run first. Phases are gated on purpose — don't skip them.
## 🔬 Pipeline overview

The analysis follows a phased evidence chain. Each phase consumes the artifacts of the previous phase and gates progression with a Go/No-Go check.

| # | Phase | Goal | Sub-skill | Artifact produced |
|:--|:------|:-----|:----------|:------------------|
| **1** | `repeat-annotation` | Soft-mask repeats (EarlGrey / RepeatModeler2 + RepeatMasker) | `repeat-annotation/` | `$RUN_DIR/assembly.fa.masked` |
| **2** | `gene-prediction-decision-tree` | Route gene prediction by evidence type and clade | `gene-prediction-decision-tree/` | routing decision recorded in the run dialogue |
| **3** | `braker3-pipeline` | Predict genes with BRAKER3 (RNA-seq and/or protein evidence) | `braker3-pipeline/` | `$RUN_DIR/braker.gff3` + `proteins.faa` |
| **4** | `galba-pipeline` | Predict genes with GALBA (protein-evidence, close relatives) | `galba-pipeline/` | `$RUN_DIR/galba.gff3` + `proteins.faa` |
| **5** | `funannotate-pipeline` | Predict genes with Funannotate (fungal clade) | `funannotate-pipeline/` | `$RUN_DIR/funannotate.gff3` + `proteins.faa` |
| **6** | `ncrna-annotation` | Annotate ncRNAs (Infernal/Rfam, tRNAscan-SE, barrnap) | `ncrna-annotation/` | combined ncRNA GFF3 |
| **7** | `functional-annotation` | Assign GO/KEGG/Pfam/EC (eggNOG-mapper + InterProScan) | `functional-annotation/` | annotated proteins + GAF + eggNOG/InterPro tables |
| **8** | `annotation-qc` | BUSCO genome vs proteome diagnostic, OMArk, sanity panel | `annotation-qc/` | `annotation-qc-report.md` |
| **9** | `runners/nextflow-runner` | Reproducible pipeline execution of Phases 1–5 | `runners/nextflow-runner/` | `results/` from `nextflow run` |

> [!TIP]
> Read each sub-skill's `SKILL.md` for the full procedure, its ask-user stop points, and its signature library (stderr pattern → cause → fix).

## 🧰 Tools

All tools are resolved from conda-forge/bioconda via the pinned `pixi.toml`.

| Tool | Version | Role |
|:-----|:--------|:-----|
| `busco` | * | qc-completeness |
| `compleasm` | * | qc-completeness |
| `earlgrey` | * | repeat-discovery-masking |
| `repeatmodeler` | * | repeat-discovery |
| `repeatmasker` | * | repeat-masking |
| `braker3` | * | gene-prediction |
| `galba` | * | gene-prediction-protein — ⚠️ not in pixi.toml (not conda-packaged); install `galba.pl` from [Gaius-Augustus/GALBA](https://github.com/Gaius-Augustus/GALBA); its deps (`augustus`, `miniprot`, `diamond`) are in the env |
| `funannotate` | * | gene-prediction-fungi |
| `infernal` | * | ncrna-cmscan |
| `trnascan-se` | * | ncrna-trna |
| `barrnap` | * | ncrna-rrna |
| `eggnog-mapper` | * | functional-orthology |
| `interproscan` | * | functional-domains |
| `helixer` | * | gene-prediction-abinitio — ⚠️ not in pixi.toml (no conda/PyPI package); install via [helixer-docker](https://github.com/gglyptodon/helixer-docker) (recommended) or manual Linux install; used only by the no-evidence route of the decision tree |
| `purge_dups` | * | haplotig-purge |
| `hisat2` | * | rnaseq-aligner |
| `samtools` | * | utility |
| `stringtie` | * | rnaseq-transcript |
| `gffread` | * | utility |
| `miniprot` | * | protein-alignment |
| `diamond` | * | protein-search |
| `augustus` | * | predictor-backend |

Tool usage is grounded in the offline `docs-corpus/` snapshots (version-matched `--help`/`man` captures, upstream repo docs, and web docs as fallback) — the agent reads these instead of guessing flags.

## 🔄 Update check

Every skill built with [skill-forge](https://github.com/cheahhl814/skill-forge) ships a self-update check that compares the deployed git SHA against this upstream repo via `git fetch` — no GitHub API call, no extra dependencies.

```bash
pixi run update-check
```

| Verdict | Exit | Meaning |
|:--------|:-----|:--------|
| `UP-TO-DATE` | 0 | Local HEAD matches origin/HEAD |
| `LOCAL-AHEAD` | 0 | Unpushed local commits; no action needed |
| `BEHIND-BY-N` | 1 | Upstream is N commits ahead → re-pull/re-sync from this repo |
| `OFFLINE` | 2 | `git fetch` failed; informational only |
| `NO-ORIGIN` | 2 | No `origin` remote configured; informational only |

## 📁 Repository layout

```text
eukaryotic-genome-annotation/
├── SKILL.md                 # Master orchestrator (router — start here)
├── README.md                # This file
├── pixi.toml                # Pinned tool environment (pixi install)
├── docs-corpus/             # Offline snapshots of upstream tool docs
│   └── busco/
│   └── compleasm/
│   └── earlgrey/
│   └── repeatmodeler/
│   └── repeatmasker/
│   └── braker3/
│   └── galba/
│   └── funannotate/
│   └── infernal/
│   └── trnascan-se/
│   └── barrnap/
│   └── eggnog-mapper/
│   └── interproscan/
│   └── helixer/
│   └── purge_dups/
│   └── hisat2/
│   └── samtools/
│   └── stringtie/
│   └── gffread/
│   └── miniprot/
│   └── diamond/
│   └── augustus/
├── repeat-annotation/            # Phase sub-skill
├── gene-prediction-decision-tree/            # Phase sub-skill
├── braker3-pipeline/            # Phase sub-skill
├── galba-pipeline/            # Phase sub-skill
├── funannotate-pipeline/            # Phase sub-skill
├── ncrna-annotation/            # Phase sub-skill
├── bin/
│   └── skill-update-check.py  # Self-update check (pixi run update-check)
└── bin/
    └── scaffold-render.py     # Reproducibility — re-render the skill from params.json
```

## 🔒 Hard guarantees

- **Filesystem evidence chain** — every phase emits the artifact the next phase consumes; the boundary between phases is the filesystem, not agent memory.
- **Docs-grounded tool usage** — tool flags come from `docs-corpus/`, never from the agent's memory.
- **Explicit stop points** — ambiguous decisions are surfaced to you as *Evidence + Recommend + Options*, not auto-picked.
- **Reproducible environments** — every tool is pinned in `pixi.toml` and resolved via pixi.

## 🚀 Nextflow runner

For production / HPC / cohort runs, the skill ships a thin Nextflow DSL2 runner at `runners/nextflow-runner/` (nf-core-style resource labels, pinned containers, trace/report/timeline enabled).

```bash
# Stub smoke test (builds the DAG without executing)
nextflow run runners/nextflow-runner/main.nf -profile test -stub-run

# Real run
nextflow run runners/nextflow-runner/main.nf --input <samplesheet.csv> --outdir <outdir> -profile docker
```

The bash recipes remain the source of truth; the runner is a thin executor.

## Provenance

Built with the [skill-forge](https://github.com/cheahhl814/skill-forge) meta-skill (v1.1.0), following the AiX-BIO skill convention (preflight → build → debug → battle-test evidence chain). Pattern adopted from:

- **BettaMt-agents** — https://github.com/cheahhl814/BettaMt-agents
- **bacterial-genome-analysis** — https://github.com/cheahhl814/bacterial-genome-analysis
- **amr-gene-screening** — https://github.com/cheahhl814/amr-gene-screening

## License

Released under the MIT License — see the license file in this repository for details.