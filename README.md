# BBB Permeability Prediction Using Spectral Topological Descriptors

This repository contains the Python implementation and the calculated spectral descriptor dataset used in the research manuscript:

> **Predicting Blood-Brain Barrier Permeability Using Eigenvalue-Based Spectral Topological Descriptors and Machine Learning**

The repository is intended to support the reproducibility of the spectral descriptor calculations presented in the manuscript.

---

# Overview

Blood–brain barrier (BBB) permeability prediction plays an important role in central nervous system (CNS) drug discovery. This repository provides the Python implementation developed to calculate fourteen eigenvalue-based spectral topological descriptors from molecular graphs. These descriptors were subsequently used for machine learning analysis in the accompanying study.

The program constructs molecular graphs from SMILES strings, computes various graph matrices, determines their eigenvalues, and calculates fourteen spectral descriptors.

---

# Repository Contents

```
BBB-Permeability-Prediction/
│
├── Descriptor_Calculation.py
│      Python implementation for spectral descriptor calculation
│
├── BBBP_14_Descriptors.xlsx
│      Calculated spectral descriptor dataset
│
└── README.md
       Repository documentation
```

---

# Input

The program accepts molecular structures represented as **SMILES** strings.

Example

```
CC(=O)Oc1ccccc1C(=O)O
```

---

# Output

For every molecular graph, the program calculates the following fourteen eigenvalue-based spectral topological descriptors.

| Symbol | Descriptor |
|---------|------------|
| EA | Adjacency Energy |
| ρA | Adjacency Spectral Radius |
| EL | Laplacian Energy |
| ρL | Laplacian Spectral Radius |
| EQ | Signless Laplacian Energy |
| ρQ | Signless Laplacian Spectral Radius |
| EES | Eccentricity Matrix Energy |
| ρES | Eccentricity Matrix Spectral Radius |
| ERS | Randic Matrix Energy |
| ρRS | Randic Matrix Spectral Radius |
| ESS | Seidel Matrix Energy |
| ρSS | Seidel Matrix Spectral Radius |
| EGAS | Graph Arithmetic–Geometric Matrix Energy |
| ρGAS | Graph Arithmetic–Geometric Matrix Spectral Radius |

---

# Dataset

The repository includes the calculated spectral descriptor dataset generated from the MoleculeNet BBBP dataset.

Each record contains

- SMILES representation
- BBB permeability label
- Fourteen spectral topological descriptors

The descriptor dataset was generated using the accompanying Python implementation.

---

# Methodology

The descriptor calculation pipeline consists of the following steps.

1. Read molecular structures from the input dataset.
2. Convert each SMILES string into a molecular graph using RDKit.
3. Construct graph matrices.
4. Compute eigenvalues of each graph matrix.
5. Calculate fourteen eigenvalue-based spectral descriptors.
6. Store the calculated descriptors in a tabular dataset for subsequent machine learning analysis.

---

# Software Requirements

The code was developed using

- Python 3
- RDKit
- NumPy
- Pandas
- tqdm

The following Python standard libraries are also used

- warnings
- urllib.request

Install the required packages using

```bash
pip install numpy pandas rdkit tqdm
```

---

# Usage

Run

```bash
python Descriptor_Calculation.py
```

The program will

- Read molecular structures
- Calculate spectral descriptors
- Generate the descriptor dataset

---

# Scientific Background

The implemented descriptors are based on spectral graph theory and are computed from eigenvalues of several graph matrices, including

- Adjacency Matrix
- Laplacian Matrix
- Signless Laplacian Matrix
- Eccentricity Matrix
- Randic Matrix
- Seidel Matrix
- Graph Arithmetic–Geometric Matrix

These descriptors provide numerical representations of molecular structures for quantitative structure–property relationship (QSPR) modelling and machine learning applications.

---

# License

This repository is made available for academic and research purposes.

Please cite the associated publication when using this code or dataset.

---

# Contact

For questions regarding the code or dataset, please contact the corresponding author through the details provided in the associated publication.
