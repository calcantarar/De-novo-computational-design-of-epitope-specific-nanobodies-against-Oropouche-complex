This document details the configuration parameters used across the computational workflow for the *de novo* nanobody design targeting the **Orthobunyavirus oropoucheense** nucleocapsid (N) protein.

---

## Defined Epitope Hotspots

| Epitope | Target Residues | Target Description |
| :--- | :---: | :--- |
| **Epitope B** | `A78`, `A79`, `A81`, `A84` | Residues 75–85 (Pan-specific & conserved) |
| **Epitope C** | `A176`, `A180`, `A182`, `A185` | Residues 171–186 (Pan-specific & conserved) |

---

## De Novo Design Configurations

### Scaffold 1: `9EZU`

<details open>
<summary><b>Click to expand/collapse Scaffold 9EZU parameters</b></summary>

```yaml
# ------------------------------------------------------------------------------
# 1. RFdiffusion Protocol
# ------------------------------------------------------------------------------
name: "orovn_9ezu"
pdb: "cluster.pdb"
hotspot: "(Defined per target epitope)"
num_designs: 8
iterations: 50
contigs: "A1-231/0 B1-23/7-13/B33-49/6-12/B58-97/12-22/B113-124"

# ------------------------------------------------------------------------------
# 2. ProteinMPNN Redesign
# ------------------------------------------------------------------------------
num_seqs: 8
mpnn_sampling_temp: 0.1
use_solubleMPNN: TRUE
rm_aa: "C" # Exclude Cysteines to prevent unwanted disulfide bridges

# ------------------------------------------------------------------------------
# 3. AlphaFold Validation
# ------------------------------------------------------------------------------
num_recycles: 3
use_multimer: TRUE
initial_guess: FALSE
```

</details>

---

### 🔶 Scaffold 2: `8QF4`

<details open>
<summary><b>Click to expand/collapse Scaffold 8QF4 parameters</b></summary>

```yaml
# ------------------------------------------------------------------------------
# 1. RFdiffusion Protocol
# ------------------------------------------------------------------------------
name: "orovn_8qf4"
pdb: "cluster.pdb"
hotspot: "(Defined per target epitope)"
num_designs: 100-500
iterations: 50
contigs: "A1-231/0 B1-25/8-14/B36-52/6-12/B61-100/12-22/B117-128"

# ------------------------------------------------------------------------------
# 2. ProteinMPNN Redesign
# ------------------------------------------------------------------------------
num_seqs: 8
mpnn_sampling_temp: 0.1
use_solubleMPNN: TRUE
rm_aa: "C" # Exclude Cysteines to prevent unwanted disulfide bridges

# ------------------------------------------------------------------------------
# 3. AlphaFold Validation
# ------------------------------------------------------------------------------
num_recycles: 3
use_multimer: TRUE
initial_guess: FALSE
```

</details>
