# Environment

Machine and software versions measured on the analysis workstation (recorded 2026-08-30).

## Machine

| Item | Value |
|---|---|
| System | macOS 26.6.2 (Build 25G83), arm64 |
| CPU | Apple M2 Pro, 12 cores |
| Memory | 32 GB unified memory |

## R

- Default R: **4.5.3 (2026-03-11)**
- R versions are managed with [rig](https://github.com/r-lib/rig), which keeps multiple R versions side by side:
  - `4.3-arm64` (R 4.3.3) — Seurat 4.4.0 / SeuratObject 4.1.4 (spatial transcriptomics analysis)
  - `4.5-arm64` (R 4.5.3) — default (bulk RNA-seq, machine learning, SHAP)

### Key R packages (R 4.5.3)

| Package | Version |
|---|---|
| Mime1 | 0.0.0.9000 (GitHub development build) |
| limma | 3.66.0 |
| ConsensusClusterPlus | 1.74.0 |
| clusterProfiler | 4.18.4 |
| fastshap | 0.1.1 |
| randomForestSRC | 3.3.2 |
| survival | 3.8.6 |
| timeROC | 0.4.1 |
| oncoPredict | 1.3.1 |
| Seurat | 5.5.0 |
| ggplot2 | 4.0.3 |
| rms | 8.1.1 |
| rmda | 1.6 |
| immunedeconv | 2.1.0 |
| scTenifoldKnk | 1.0.3 |
| scatterplot3d | 0.3.45 |
| AnnotationDbi | 1.72.0 |

## Molecular dynamics / cheminformatics

Dedicated mamba environment (`md-clean`):

| Tool | Version |
|---|---|
| GROMACS | 2025.4 (conda-forge build) |
| AmberTools | 24.8 |
| ACPYPE | 2023.10.27 |
| Open Babel | 3.1.1 |
| Python | 3.12.13 |

## Python environments (miniconda)

| Environment | Purpose | Key versions |
|---|---|---|
| graphban | GraphBAN drug-screening model | Python 3.11.15; torch 2.3.1; torch-geometric 2.6.1; pytorch-lightning 2.6.1; RDKit 2024.3.6 |
| graphban-clean | Cleaned variant of graphban | Python 3.11.15; torch 2.3.1; torch-geometric 2.6.1; RDKit 2026.3.1 |
| admet-ai | ADMET property prediction | Python 3.11.15; admet_ai 2.0.1 |
| md-2025.4 | MD trajectory analysis | Python 3.12.14; MDAnalysis 2.10.0; GROMACS 2025.4 |

## Project conventions

- Random seed: **2026** for all stochastic steps
- Colors: tumor `#E64B35`, normal `#4DBBD5`
- SHAP: `fastshap::explain()` with `nsim = 200`, 4 parallel cores
- Spatial transcriptomics analysis uses Seurat v4 (R 4.3.3); the main analysis uses Seurat v5 (R 4.5.3)
