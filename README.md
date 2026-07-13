# BBB Permeability Prediction

This repository contains the Python implementation used to calculate the fourteen eigenvalue-based spectral topological descriptors employed in the manuscript:

**Predicting Blood-Brain Barrier Permeability Using Eigenvalue-Based Spectral Topological Descriptors and Machine Learning**

## Overview

The program calculates fourteen eigenvalue-based spectral descriptors from molecular graphs constructed from SMILES strings.

## Calculated Spectral Descriptors

- Adjacency Energy (EA)
- Adjacency Spectral Radius (ρA)
- Laplacian Energy (EL)
- Laplacian Spectral Radius (ρL)
- Signless Laplacian Energy (EQ)
- Signless Laplacian Spectral Radius (ρQ)
- Eccentricity Matrix Energy (EES)
- Eccentricity Matrix Spectral Radius (ρES)
- Randic Matrix Energy (ERS)
- Randic Matrix Spectral Radius (ρRS)
- Seidel Matrix Energy (ESS)
- Seidel Matrix Spectral Radius (ρSS)
- Graph Arithmetic–Geometric Matrix Energy (EGAS)
- Graph Arithmetic–Geometric Matrix Spectral Radius (ρGAS)

## Programming Language

- Python 3

## Main Library Requirements

- RDKit
- NumPy
- NetworkX
- Pandas

This repository is provided to support the reproducibility of the spectral descriptor calculations presented in the accompanying manuscript.
