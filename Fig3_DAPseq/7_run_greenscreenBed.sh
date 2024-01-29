#!/bin/bash
#SBATCH -p nbi-medium
#SBATCH -t 0-04:00
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/peaks/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/peaks/slurm_output/%x.%N.%j.err

cd /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks

#bedtools 2.29.2
source package 6394519c-541f-479c-b064-dd0b912eac04

echo Packages loaded. >&2

#qval
#merge_distance
#n_distinct inputs

srun bash scripts/generate_3input_greenscreenBed.sh 5 20000 2
srun bash scripts/generate_3input_greenscreenBed.sh 5 30000 2
srun bash scripts/generate_3input_greenscreenBed.sh 5 50000 2
srun bash scripts/generate_3input_greenscreenBed.sh 5 100000 2

srun bash scripts/generate_3input_greenscreenBed.sh 10 20000 2
srun bash scripts/generate_3input_greenscreenBed.sh 10 30000 2
srun bash scripts/generate_3input_greenscreenBed.sh 10 50000 2
srun bash scripts/generate_3input_greenscreenBed.sh 10 100000 2

srun bash scripts/generate_3input_greenscreenBed.sh 2 20000 2
srun bash scripts/generate_3input_greenscreenBed.sh 2 30000 2
srun bash scripts/generate_3input_greenscreenBed.sh 2 50000 2
srun bash scripts/generate_3input_greenscreenBed.sh 2 100000 2


