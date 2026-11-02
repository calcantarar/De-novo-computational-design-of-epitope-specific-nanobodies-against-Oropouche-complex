# *De novo* Computational Design of Epitope-Specific Nanobodies Against the *Orthobunyavirus oropoucheense* Nucleocapsid for Point-of-Care Diagnostics

\---

## Abstract

### Background

The 2023–2025 Oropouche outbreak, with over 23,000 cases across the American continent, exposed critical public health gaps. Although RT-qPCR remains the diagnostic gold standard, it requires specialized infrastructure inaccessible in endemic rural areas. Coupled with reinfection rates approaching 60%, these gaps underscore the urgent need for low-cost point-of-care (POC) alternatives such as lateral flow assays. Traditional nanobody discovery relies on animal immunization, which is slow and unable to target predefined epitopes, and no epitope-directed POC reagents have been reported for the *Orthobunyavirus oropoucheense* species complex.

### 

### Methodology

We developed a multi-stage *in silico* workflow to design epitope-specific nanobodies against the nucleocapsid (N) protein of this species complex:

1. **Sequence \& Structure Curation:** Nucleocapsid sequences spanning OROV and closely related strains were curated via **MMseqs2** and modeled with **AlphaFold2 / ColabFold**.
2. **Epitope Selection:** Linear and discontinuous epitopes were predicted using **DTU Bioinformatic Services** and filtered by strict sequence conservation and minimal RNA-binding interference.
3. *De novo* **Design \& Docking:** Validated VHH scaffolds (PDB: `8QF4`, `9EZU`) underwent forced docking with **HADDOCK v2.4**, followed by *de novo* CDR redesign using **RFdiffusion**, **ProteinMPNN**, and **AlphaFold-Multimer**.
4. **Filtering \& MD Validation:** Candidates were progressively filtered through structural metrics, **PDBePISA** interfacial ranking, and 100 ns all-atom molecular dynamics (**GROMACS 2023**, AMBER99SB-ILDN, TIP3P).

### 

### Results \& Key Findings

* **Identified Epitopes:** Two pan-specific epitopes were selected:

  * **Epitope B** (residues 75–85)
  * **Epitope C** (residues 171–186)
  * *Both showed high conservation across the Oropouche complex and zero cross-reactivity against Zika, Dengue, Chikungunya, and Mayaro viruses (BLASTp, E-value > 1.0).*
* **Design Selection:** Out of 96 generated designs, 14 passed initial structural filters.

  * **Epitope C Candidates:** Met strict metrics ($\\text{pLDDT} \\ge 0.75$, $\\text{pTM} \\ge 0.70$, $\\text{PAE} \\le 10\\text{ \\AA}$, $\\text{RMSD} \\le 2.5\\text{ \\AA}$); top complex achieved $\\text{pLDDT} = 0.863$, $\\text{pTM} = 0.793$, and interface $\\text{RMSD} = 1.03\\text{ \\AA}$.
  * **Epitope B Candidates:** Retained under relaxed criteria due to higher flexibility (best: $\\text{pLDDT} = 0.78$, $\\text{pTM} = 0.57$, $\\text{RMSD} = 4.2\\text{ \\AA}$).
* **Interface Interactions (PDBePISA):**

  * **Epitope B:** 4 H-bonds and 1 salt bridge.
  * **Epitope C:** 3 H-bonds and 3 salt bridges.
* **100 ns MD Simulations:** Sustained an average of **7.66** (Epitope B) and **7.67** (Epitope C) intermolecular H-bonds per frame. Despite backbone RMSD reaching 7.98 Å and 8.32 Å due to conformational sampling, the persistent interfacial H-bond network indicates the binding interface remained continuously engaged.

### Conclusions

This workflow provides a generalizable, epitope-directed framework for rapid, low-cost nanobody design against emerging arboviruses, supporting the accelerated development of species-specific POC lateral flow assays for resource-limited settings in Latin America.

\---

**Keywords:** *De novo* protein design · Nanobodies · *Orthobunyavirus oropoucheense* · Point-of-care diagnostics · Molecular dynamics

