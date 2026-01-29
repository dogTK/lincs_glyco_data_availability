# LINCS L1000 Glycogene Analysis

Analysis code and results for the study of glycosylation-related gene expression patterns in the LINCS L1000 chemical perturbation dataset.

## Overview

This repository contains the analysis pipeline and results for characterizing drug-induced transcriptional responses of glycosylation-related genes (glycogenes) using the LINCS L1000 dataset.

### Key Features

- **Data Extraction**: Scripts to extract glycogene expression data from LINCS L1000 Chemical Perturbations (2021) dataset
- **Rank-Based Analysis**: Pure rank aggregation method for identifying consistently responsive glycogenes
- **WGCNA Analysis**: Weighted gene co-expression network analysis to identify glycogene modules
- **ATC Classification**: Integration with WHO ATC drug classification for therapeutic class enrichment analysis
- **PHH Similarity Analysis**: Comparison with primary human hepatocyte expression patterns

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
│   │   ├── 09c_atc_pathway_gsea.ipynb
│   │   └── 10_phh_rank_similarity_analysis.ipynb
│   ├── engineering/        # Data preparation notebooks
│   │   ├── 00_create_glyco_genes_wide_table.ipynb
│   │   ├── 00b_create_atc_master_table.ipynb
│   │   ├── 01_create_ctl_rnaseq_glyco_table.ipynb
│   │   └── 11_compound_name_mapping.ipynb
│   ├── GlycoEnzOnto/       # Glycosylation enzyme ontology
│   └── results/            # Analysis outputs
├── src/                    # Python modules
│   ├── compound_name_mapping.py
│   ├── task2_phh_rank_similarity.py
│   └── utils/
└── environment.yml         # Conda environment specification
```

## Data Sources

- **LINCS L1000**: [LINCS Data Portal](https://lincsportal.ccs.miami.edu/)
- **GlycoEnzOnto**: Glycosylation enzyme ontology for pathway definitions
- **ChEMBL**: ATC classification data for drug annotation

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
2. Configure Snowflake credentials (for data extraction from private database)
3. Run notebooks in order:
   - `engineering/` notebooks for data preparation
   - `analytics/` notebooks for analysis

Note: Some notebooks require access to a Snowflake database with pre-loaded LINCS data.

## Citation

If you use this code, please cite:

> [Paper citation to be added upon publication]

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contact

For questions or issues, please open a GitHub issue.
