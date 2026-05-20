# I-indigotica-fruit-metabolomics
R scripts used for OPLS-DA, VIP extraction, and DAM analysis in Isatis indigotica fruit metabolomics.
# R scripts for *Isatis indigotica* fruit metabolomics analysis

This repository contains R scripts used for OPLS-DA, VIP extraction, and differentially accumulated metabolite (DAM) statistics in the revised manuscript:

**Stage-resolved multi-omics links flavonoid accumulation to antioxidant capacity in *Isatis indigotica* fruit**

## Purpose

The repository was created to improve computational transparency and reproducibility of the metabolomics analysis described in the manuscript.

The main script performs:

- data import
- extraction of sample columns
- log2 transformation
- pairwise OPLS-DA using the `ropls` package
- 200-permutation testing
- VIP extraction
- OPLS-DA score plot generation
- fold-change calculation
- two-sided Student's *t*-test
- Benjamini-Hochberg FDR correction
- DAM identification using VIP > 1, |log2FC| >= 1, and FDR < 0.05

## Repository structure

```text
I-indigotica-multiomics-scripts/
├── README.md
├── scripts/
│   └── 01_OPLS_DA_ropls_and_DAM.R
├── example/
│   └── example_input_format.csv
├── sessionInfo.txt
├── .gitignore
└── LICENSE
