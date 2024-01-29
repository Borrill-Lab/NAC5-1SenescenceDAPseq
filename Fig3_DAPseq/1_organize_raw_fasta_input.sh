#!/bin.bash

directory="/jic/scratch/groups/Philippa-Borrill/Catherine/raw_data/X204SC23033228-Z01-F002/X204SC23033228-Z01-F002/01.RawData"

outdir="/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/raw"
mkdir -p ${outdir}

#DAP_53
for pair in {1..2}; do
	x="DAP_53"
	#Add two files where there are two
	gzip_f1="${directory}/${x}/DAP_53_EKDL230009916-1A_HGW5VDSX7_L2_${pair}.fq.gz"
	gzip_f2="${directory}/${x}/DAP_53_EKDL230009916-1A_HKTC2DSX7_L4_${pair}.fq.gz"
	gzip_f="${outdir}/input${x}_${pair}.fastq.gz"

	echo ${gzip_f1} ${gzip_f}
	
	cat ${gzip_f1} ${gzip_f2} \
	> ${gzip_f}
done

#DAP_54
for pair in {1..2}; do
        x="DAP_54"
        #Add two files where there are two
        gzip_f1="${directory}/${x}/DAP_54_EKDL230009925-1A_HFJJGDSX7_L1_${pair}.fq.gz"
        gzip_f2="${directory}/${x}/DAP_54_EKDL230009925-1A_HKTC2DSX7_L2_${pair}.fq.gz"
        gzip_f="${outdir}/input${x}_${pair}.fastq.gz"

        echo ${gzip_f1} ${gzip_f}

        cat ${gzip_f1} ${gzip_f2} \
        > ${gzip_f}
done

directory="/jic/scratch/groups/Philippa-Borrill/Catherine/raw_data/X204SC23033228-Z01-F001/X204SC23033228-Z01-F001_01/01.RawData"

#DAP_52
for pair in {1..2}; do
        x="DAP_52"
        #Add two files where there are two
        gzip_f1="${directory}/${x}/DAP_52_EKDL230009915-1A_HGW5VDSX7_L2_${pair}.fq.gz"
        gzip_f="${outdir}/input${x}_${pair}.fastq.gz"

        echo ${gzip_f1} ${gzip_f}

        cat ${gzip_f1} \
        > ${gzip_f}
done

