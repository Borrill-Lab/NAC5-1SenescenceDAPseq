#!/bin/bash
#SBATCH -p nbi-short
#SBATCH -t 0-01:00
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/motifs/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/motifs/slurm_output/%x.%N.%j.err

cd /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks

#bedtools 2.29.2
source package 6394519c-541f-479c-b064-dd0b912eac04

>&2 echo Packages loaded.

bedtools getfasta -fi /jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta  \
-bed peaks/NAMA1/NAMA1_peaks.narrowPeak -fo motifs/NAMA1_peaks.fasta
