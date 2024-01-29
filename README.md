# NAC5-1SenescenceDAPseq
This repository contains all scripts used to generate figures in:

Evans et al. Wheat NAC transcription factor _NAC5-1_ is a positive regulator of senescence

## Contents of this repository

`Fig1_NAC5_TILLING` contains R scripts to calculate senescence metrics and plot Fig 1 and Fig S3.

`Fig2_NAC5_transgenics` contains R scripts to calculate senescence metrics, analyse qPCR data and plot Fig 2, Fig S4 and Fig S5.

`Fig3_DAPseq` contains a bash, perl and R pipeline to process DAP-seq data, plot Fig 3 and generate Supplementary Datasets.

`FigS2` contains a script to plot expression data in Fig S2.

Within each folder, scripts are numbered in the order in which they were run. Scripts without numbers are sourced indirectly.

## DAP-seq data availability
Raw data for the DAP-seq can be obtained through BioProject ID PRJEB72016 on the European Nucleotide Archive.
Within the Fig3_DAPseq folder, subfolder 'metadata' contains metadata files such as sample lists input into scripts.
The wheat reference genome Chinese Spring IWGSC Refseq v1.1 and high confidence gene annotations were obtained from https://wheat-urgi.versailles.inra.fr/Seq-Repository/

## Acknowledgement
The DAP-seq analysis pathway in this repository is inspired by (and script "generate_3input_greenscreenBed.sh" directly adapted from) the following repo, to which these authors are grateful:
Klasfeld S, Roulé T, Wagner D (2022) Greenscreen: A simple method to remove artifactual signals and enrich for true peaks in genomic datasets including ChIP-seq data. Plant Cell. doi:10.1093/plcell/koac282
https://github.com/sklasfeld/GreenscreenProject
