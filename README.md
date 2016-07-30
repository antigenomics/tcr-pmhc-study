## Meta-analysis of T-cell receptor:peptide:MHC (TCR:pMHC) complex structural data available in PDB

> **DISCLAIMER** You are free to use the processed and annotated TCR:pMHC tables in this repo for your research. However, the authors request priority in publishing a meta-analysis of the dataset aggregated and stored here. This message will be removed once this study is officially published. Thank you for understanding.

This repository contains a set of scripts that automate downloading, pre-processing and annotation of TCR:pMHC structural data with the following functionality:

**scripts in ``annotation/src/`` folder**

- Extracts PDB entry metadata and verifies that a given PDB entry represents a valid TCR:pMHC complex
- Annotates each molecule in the complex: determines MHC class and allele, peptide molecule and TCR chains
- Performs V/D/J mapping for TCR chains, partitions TCR chain into CDR/FR regions

**scripts in ``structure/src/`` folder**

- Computes pairwise amino acid distances and point energies for TCR and peptide residues

To run the pipeline execute the ``run.sh`` script. It will proceed with the list of PDB ids from ``result/extended_pdb_ids.txt``. Some parts of the script are rather time-consuming, especially structural data annotation (downloading PDB files and running GROMACS). The results will be stored in ``result/`` folder:

- ``final.annotations.txt`` contains the list of PDB entries that passed filtering and their annotation
- ``structure.txt`` contains annotated data on the amino acid level, with pairwise residue distances and interaction energies

### Pre-requisites

The pipeline is written in [Groovy](http://www.groovy-lang.org) and [Python](https://www.continuum.io/downloads) (written in ``3.5`` but should run under ``2.7``) and requires both to run. Python scripts require ``pandas`` and ``BioPython`` packages to run.

Three third-party software tools that are required:

- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). Ensure that ``blastp`` and ``makeblastdb`` are in your ``$PATH``.
- [IgBlast](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone), strictly the ``1.4.0`` version. Ensure that ``igblastp`` is in your ``$PATH``.
- [GROMACS](http://www.gromacs.org) for computing interaction energies. Ensure that ``gmx`` is in your ``$PATH``.

