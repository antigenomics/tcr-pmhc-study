# Pairwise contact models (PCMs)

This folder contains various TCR-pMHC PCMs. The models are separated into local and global:

* **Local** models only rely on residue index and CDR/antigen length to estimate pairwise residue distance and select residue pairs eligible for contact
* **Global** models use full CDR and antigen sequences to predict residue positions in 3D, rotate CDR against antigen and compute distances

Currently, all models rely on a generalized linear model (**GLM**) to predict contact probability given pairwise residue distance and types of the residues