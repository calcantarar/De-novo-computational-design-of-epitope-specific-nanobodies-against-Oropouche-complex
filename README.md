# De novo computational design of epitope-specific nanobodies against the _Orthobunyavirus oropoucheense_ nucleocapsid for point-of-care diagnostics
Alcántara-Ríos, C. ; Hernandez, A. ; Sandoval, G.

## Overview
This repository contains the computational pipeline, structural models, analysis scripts, and simulation protocols used for the *de novo* design and dynamic validation of pan-specific nanobodies (VHHs) targeting the nucleocapsid (N) protein of the *Orthobunyavirus oropoucheense* species complex.

## Repository Structure
```text
├── data/
│   ├── raw_sequences/        # MMSeqs2 curated sequence alignments
│   └── pdb_models/           # Initial and refined PDB structures (Epitope B & C complex designs)
├── protocols/
│   ├── rfdiffusion/          # Config files for CDR redesign
│   └── gromacs/              # .mdp files for energy minimization, NVT, NPT, and 100 ns production
├── figures/
│   └── MD_figures.py         # Python script for RMSD, RMSF, and H-bond extraction
├── docs/                     # Poster abstract, figures, and supplemental tables
└── README.md
