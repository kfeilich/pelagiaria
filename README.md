# pelagiaria

This directory contains code for the paper [A phylogenomic framework for pelagiarian fishes (Acanthomorpha: Percomorpha) highlights mosaic radiation in the open ocean](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2019.1502). If you want to run this code, I recommend downloading this code and the data directory from [Data Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.sq067ng).

## r-code : The directory containing the R scripts used to conduct the analyses. These should be sourced in the order in which they are numbered.
- 00_start-up.R : Imports required libraries with some supplemental functions from github (REF).
- 01_functions.R : Contains custom functions written for downstream analyses in this paper (including function documentation).
- 02_DTT_analysis.R : Runs the modified MDI analysis and the rank-envelope tests of disparity as defined by Murrell, D. J. 2018.
- 03_plot_DTT_figures.R : Calculates the MDI from the output of the DTT analysis, and plots relevant results.
- 04_PelagiariaPCA.R : Imports and cleans landmarks from TPS file, calculates taxon mean landmarks from species with multiple representatives, conducts generalized procrustes analysis and principal components analysis on non-hairtail taxa. 
- 04b_Pelagiaria_w_hairtails_PCA.R : As above, but conducts analysis on a subset of landmarks that are suitable for use with hairtail specimens, in addition to the rest. 
- 05_variance_analysis.R : measures the proportion of morphological variance present in the time window of interest relative to the morphological variance present in modern taxa. 
