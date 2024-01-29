#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 0-08:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=end,fail
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/slurm_output/%x.%N.%j.err

cd /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks

#bowtie2 2.4.1
source package 29a74b59-88fc-4453-a30b-1310b34910b9

#create directory
mkdir -p ./meta/WheatGenome/bowtie2_genome_dir

#Index genome for Bowtie2 (run once before bowtie2 can be used)
bowtie2-build --threads 8 --large-index \
/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta \
meta/WheatGenome/bowtie2_genome_dir/iwgsc_refseq_1.1
