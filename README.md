# LINCS L1000 Glycogene Analysis

Analysis code and results for the study of glycosylation-related gene expression patterns in the LINCS L1000 chemical perturbation dataset.

## Overview

This repository contains the analysis pipeline and results for characterizing drug-induced transcriptional responses of glycosylation-related genes (glycogenes) using the LINCS L1000 dataset.

### Key Features

- **Data Extraction**: Scripts to extract glycogene expression data from LINCS L1000 Chemical Perturbations (2021) dataset
- **Rank-Based Analysis**: Pure rank aggregation method for identifying consistently responsive glycogenes
- **WGCNA Analysis**: Weighted gene co-expression network analysis to identify glycogene modules
- **ATC Classification**: Integration with WHO ATC drug classification for therapeutic class enrichment analysis
- **Mechanism-of-action analysis**: Module and pathway-program responses of Drug Repurposing Hub mechanism-of-action classes, including HDAC inhibitors and compounds with an HDAC inhibitor-like response

## Repository Structure

```
.
├── notebooks/
│   ├── analytics/          # Main analysis notebooks
│   │   ├── 01_data_extraction_and_qc.ipynb
│   │   ├── 02_calculate_scores_rank_pure.ipynb
│   │   ├── 07_wgcna_analysis_v2.ipynb
│   │   ├── 08_atc_module_eigengene_analysis.ipynb
│   │   ├── 09_drug_repurposing_hub_enrichment.ipynb
│   │   └── 09c_atc_pathway_gsea.ipynb
│   ├── engineering/        # Data preparation notebooks
│   │   ├── 00_create_glyco_genes_wide_table.ipynb
│   │   ├── 00b_create_atc_master_table.ipynb
│   │   ├── 01_create_ctl_rnaseq_glyco_table.ipynb
│   │   └── 11_compound_name_mapping.ipynb
│   └── GlycoEnzOnto/       # Glycosylation enzyme ontology
├── scripts/                # Analysis scripts for the manuscript figures and tables
│   ├── 20_moa_glyco_programs/       # MOA class x module / program comparisons, Figures 2 and 3A,B,D
│   ├── 22_compose_main_figures/     # assembles Figures 1-3
│   ├── 23_hdaci_cross_cell_line/    # HDAC inhibitors in nine further cell lines, Figure 3C
│   ├── 24_program_branch_analysis/  # nucleotide sugar / sialylation branches and genes
│   ├── 25_supplementary_package/    # builds the supplementary PDF and Excel files
│   ├── 26_microtubule_agents/       # microtubule-targeting agents (not used in the manuscript)
│   └── 27_hdaci_like_compounds/     # compounds with an HDAC inhibitor-like response, Figure 3E,F
├── src/                    # Python modules
│   ├── compound_name_mapping.py
│   └── utils/
└── environment.yml         # Conda environment specification
```

## Revision scripts (2026-11)

Scripts used for the revised manuscript. They expect the project layout of the analysis repository
(`results/02_wgcna/`, `results/04_phh_similarity/`, `results/12_phh_drug_class_enrichment/`, `data/external/`)
produced by the notebooks above.

| Script | Purpose |
|---|---|
| `scripts/21_phh_axis_contributing_genes/run_phh_contributing_genes.py` | Genes driving similarity to the PHH reference signature and their GlycoEnzOnto enrichment |
| `scripts/20_hdac_phh_repurposing_hub/run_hub_phh_analysis.py` | Drug Repurposing Hub annotation, ranked enrichment of mechanism-of-action classes along the PHH similarity ranking (10,000 permutations), HDAC inhibitor robustness (perturbation-strength-matched control, per-condition similarity, module and gene profiles), and program-level shifts for all classes |
| `scripts/23_hdaci_cross_cell_line/run_cross_cell_line.py` | HDAC inhibitor program-level shifts in eight additional LINCS cell lines |
| `scripts/plot_integrated_atc_phh_module_complexheatmap.R` | ATC class x module GSEA heatmap with module eigengene distributions along the PHH axis |

`run_cross_cell_line.py` queries Snowflake with key-pair authentication; set `SNOWFLAKE_USER`, `SNOWFLAKE_ACCOUNT`,
`SNOWFLAKE_PRIVATE_KEY_PATH`, and optionally `SNOWFLAKE_WAREHOUSE` and `SNOWFLAKE_ROLE`.

## Data Sources

- **LINCS L1000**: [LINCS Data Portal](https://lincsportal.ccs.miami.edu/)
- **GlycoEnzOnto**: Glycosylation enzyme ontology for pathway definitions
- **ChEMBL**: ATC classification data for drug annotation
- **Drug Repurposing Hub** (version 2020-03-24): mechanism-of-action and target annotations, https://repo-hub.broadinstitute.org/repurposing

## Requirements

```bash
conda env create -f environment.yml
conda activate l1000-glyco
```

Key dependencies:
- Python 3.11
- pandas, numpy, scipy
- matplotlib, seaborn
- scikit-learn
- cmapPy (for GCTX file parsing)

## Usage

1. Set up the conda environment
2. Configure Snowflake credentials (see below)
3. Run notebooks in order:
   - `engineering/` notebooks for data preparation
   - `analytics/` notebooks for analysis

### Snowflake Configuration

The notebooks require access to a Snowflake database with pre-loaded LINCS data. Set the following environment variables:

```bash
export SNOWFLAKE_USER="your_username"
export SNOWFLAKE_ACCOUNT="your_account_identifier"
export SNOWFLAKE_PRIVATE_KEY_PATH="~/.ssh/your_private_key.pem"
```

The code uses key-pair authentication. See [Snowflake documentation](https://docs.snowflake.com/en/user-guide/key-pair-auth) for setup instructions.

## Citation

If you use this code, please cite:

> [Paper citation to be added upon publication]

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contact

For questions or issues, please open a GitHub issue.
