[![DOI Code](https://zenodo.org/badge/DOI/10.5281/zenodo.23994366.svg)](https://doi.org/10.5281/zenodo.23994366)
[![DOI Dataset](https://img.shields.io/badge/Zenodo-MD_Dataset-blue.svg)](https://doi.org/10.5281/zenodo.23965995)

# De novo computational design of epitope-specific nanobodies against the _Orthobunyavirus oropoucheense_ nucleocapsid for point-of-care diagnostics
Alcántara-Ríos, C. ; Hernandez, A. ; Sandoval, G.

## Overview
This repository contains the computational pipeline, structural models, analysis scripts, and simulation protocols used for the *de novo* design and dynamic validation of pan-specific nanobodies (VHHs) targeting the nucleocapsid (N) protein of the *Orthobunyavirus oropoucheense* species complex.

## Repository Structure
```text
├── data/
│   ├── fasta_sequences/            # Sequences retrieved from NCBI virus and MMSeqs2 curated sequence alignments
│   └── pdb_models/                 # Initial and refined PDB structures
│   └── restraint-driven_docking/   # Selected results from HADDOCK
│   └── epitope_prediction/         # Merged results from DTU services, IEDB and conservation metrics of OROV N
│   └── directed_mutagenesis/       # Filtered and selected PDB structures
├── protocols/
│   ├── rfdiffusion/          # Configuration for CDR redesign in Difussion.py
│   └── gromacs/              # .mdp files for energy minimization, NVT, NPT, and 100 ns production
├── figures/
│   └── MD_figures.py         # Python script for RMSD, RMSF, and H-bond extraction
├── docs/                     # Poster, abstract and supplemental tables
└── README.md
```

## Cite and data availability

* **Code and computational pipeline (v1.0.0):** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23994366.svg)](https://doi.org/10.5281/zenodo.23994366)
  * Includes filter scripts, metric analysis, raw data and tables, processed `.xvg` files for the figure generation (`MD_figures.py`).

* **Molecular dynamics raw dataset (100 ns):** [![DOI](https://img.shields.io/badge/Zenodo-Dataset-blue.svg)](https://doi.org/10.5281/zenodo.23965995)
  * Contains the GROMACS complete trajectories (`.xtc`), topologies (`.top`/`.itp`) and initial structures (`.gro`/`.tpr`) for the VHH-OROVN complexes against Epitope B and C.
