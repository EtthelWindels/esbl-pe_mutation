# Estimating the within-host mutation rate of ESBL-PE using phylodynamics
This repository contains the code associated with the phylodynamic analyses performed in Windels, Aguilar-Bultet et al. Phylodynamic estimation of the within-host mutation rate of extended-spectrum beta-lactamase-producing Enterobacterales (2025, bioRXiv)

These analyses aimed at estimating the within-patient mutation rates of ESBL-producing *E. coli* and *K. pneumoniae* species complex using sequence data from longitudinally collected rectal swabs.
The raw sequencing data used in this study are available in the NCBI database under the BioProject number [PRJNA910977](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA910977).

The repository is organized as follows:
- The folder `analyses/` contains the XML files used for the phylodynamic analyses in BEAST2.
- The folder `alignments/` contains the alignments used as input for BEAST2.
- The folder `scripts/` contains R script used for post-processing of the BEAST2 output.
