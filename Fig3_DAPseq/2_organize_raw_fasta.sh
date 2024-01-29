#!/bin/bash
#SBATCH -p nbi-short
#SBATCH -t 0-01:00
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=Catherine.Evans@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/slurm_output/%x.%N.%j.err

dir1="/jic/scratch/groups/Philippa-Borrill/Catherine/raw_data/X204SC23033228-Z01-F001/X204SC23033228-Z01-F001_01/01.RawData"
dir2="/jic/scratch/groups/Philippa-Borrill/Catherine/raw_data/X204SC23033228-Z01-F001/X204SC23033228-Z01-F001_02/01.RawData"
dir3="/jic/scratch/groups/Philippa-Borrill/Catherine/raw_data/X204SC23033228-Z01-F002/X204SC23033228-Z01-F002/01.RawData"

outdir="/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/raw"
mkdir -p ${outdir}

#dir1
input_list=("DAP_5" "DAP_6" "DAP_7" "DAP_47")

for x in "${input_list[@]}"; do
	cd -e ${dir1}/${x}
	read1=$(find . -name "*1.fq.gz")
	read2=$(find . -name "*2.fq.gz")
	gzip_read1="${outdir}/${x}_1.fastq.gz"
	gzip_read2="${outdir}/${x}_2.fastq.gz"
#	>&2 echo read names for ${x} are ${read1} ${read2} 
#	>&2 echo to be saved as ${gzip_read1} ${gzip_read2}
#	cat ${read1} > ${gzip_read1}
#	cat ${read2} > ${gzip_read2}
done

#dir2
input_list=("DAP_1" "DAP_2" "DAP_3" "DAP_4" "DAP_12" \
        "DAP_13" "DAP_14" "DAP_17" "DAP_18" "DAP_20" "DAP_23" \
        "DAP_24" "DAP_25" "DAP_28" "DAP_29" "DAP_30" \
	"DAP_31" "DAP_37" "DAP_38" "DAP_39" "DAP_40" \
	"DAP_44" "DAP_45" "DAP_46" )

for x in "${input_list[@]}"; do
        cd -e ${dir2}/${x}
        read1=$(find . -name "*1.fq.gz")
        read2=$(find . -name "*2.fq.gz")
        gzip_read1="${outdir}/${x}_1.fastq.gz"
        gzip_read2="${outdir}/${x}_2.fastq.gz"
#        >&2 echo read names for ${x} are ${read1} ${read2}
#        >&2 echo to be saved as ${gzip_read1} ${gzip_read2}
#        cat ${read1} > ${gzip_read1}
#        cat ${read2} > ${gzip_read2}
done

#dir3
input_list=("DAP_10" "DAP_15" "DAP_16" "DAP_19" \
	"DAP_22" "DAP_32" "DAP_33" "DAP_34" "DAP_35" \
	"DAP_36" )

for x in "${input_list[@]}"; do
        cd -e ${dir3}/${x}
        read1=$(find . -name "*1.fq.gz" | sort)
        read2=$(find . -name "*2.fq.gz" | sort)
        gzip_read1="${outdir}/${x}_1.fastq.gz"
        gzip_read2="${outdir}/${x}_2.fastq.gz"
        >&2 echo read names for ${x} are ${read1} ${read2}
        >&2 echo to be saved as ${gzip_read1} ${gzip_read2}
        cat ${read1} > ${gzip_read1}
        cat ${read2} > ${gzip_read2}
done

