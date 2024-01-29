#!/usr/bin/perl -w

# 
#
# Aim of script is to run macs2 on DAP-seq data
# Pool samples from 3 technical replicates (without downsampling)

#### paths and references:
my $path = '/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/bowtie2_index';
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/mapped";
my $ref = "$path/iwgsc_refseq_1.1";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/peaks";

### list containing path to samples starting from $read_path_triticum
#####
# Your list of samples should look like this:
# Sample	treat1	treat2	treat3	control1	control2	control3

# Make sure the columns are tab separated


my $input_list_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/mapped";
my $input_for_macs2 = "sample_list_pooled.txt";

### greenscreen should be a tab separated interval file of greenscreen regions to mask
# chr1a_part1	123456000	123789000

my $q = 2;
my $greenscreen_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/peaks/qval$q";
my $greenscreen = "gs_merge50000bp_call2_inputs.txt";

### chromsizes should be a tab separated file with chromosome names and sizes
#chr1A_part1      471304005

my $chromsizes = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/scripts/chromosome_lengths.txt";

#######OUTPUT DIRECTORY, TMP DIRECTORY, SLURM_OUTPUT DIRECTORY MUST BE CREATED BEFORE RUNNING THE SCRIPT


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$input_list_dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_macs2") || die "couldn't open the input file $input_for_macs2!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
print "\nmy line was: $line\n";
  
#print "\nmy array: @array\n";
print "\nsample: @array[0]\n";
  
my $sample = $array[0];
my $t1 = $array[1];
my $t2 = $array[2];
my $t3 = $array[3];

my $c1 = $array[4];
my $c2 = $array[5];
my $c3 = $array[6];
  
chdir("$read_path_triticum") or die "couldn't move to specific read directory $read_path_triticum";
  
  
my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel macs2 tasks
#
#SBATCH -p nbi-short
#SBATCH -t 0-01:00
#SBATCH -c 1
#SBATCH --mem=60000
#SBATCH -J macs2
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/peaks/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/peaks/slurm_output/%x.%N.%j.err
SLURM
  
 my $tmp_file = "$output_dir/tmp/macs2.$sample";
  
  
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";
  
  print SLURM "set -e\n";
  print SLURM "mkdir -p $output_dir/$sample\n";
### this part all needs editing!!!! ###
    #MACS 2.2.7.1
    print SLURM "source package 0aea19f1-78a3-40ee-87d3-49ab33c04ff5\n";
  
    print SLURM "macs2 callpeak --bdg -t $t1 $t2 $t3 -c $c1 $c2 $c3 -f BAMPE --keep-dup auto -g 12365172330 -n $sample --outdir $output_dir/$sample"."\n";
    print SLURM ">&2 echo peaks called\n";

    # remove all peaks that do not have an
    # average base pair q-value <=10^(-${q})
    print SLURM "awk -F\"\\t\" -v q=${q}"." 'BEGIN{OFS=\"\\t\"} \$9>=q {print}' \\\n $output_dir/$sample/$sample"."_peaks.narrowPeak > $output_dir/$sample/$sample"."noMask_qval$q"."_peaks.narrowPeak\n";
  
    # remove all peaks that overlap greenscreen
    #bedtools 2.29.2
    print SLURM "source package 6394519c-541f-479c-b064-dd0b912eac04\n";
    print SLURM "bedtools intersect -v -wa \\\n -a $output_dir/$sample/$sample"."noMask_qval$q"."_peaks.narrowPeak -b $greenscreen_dir/$greenscreen > \\\n $output_dir/$sample/$sample"."gsMask_qval$q"."_peaks.narrowPeak\n";

    #convert BedGraph to BigWig
    #KentTools 1.0
    print SLURM "source package 71de5a7a-135b-417a-8de1-ede16dc52660\n";

    print SLURM "bedGraphToBigWig $output_dir/$sample/$sample"."_control_lambda.bdg $chromsizes $output_dir/$sample/$sample"."_control_lambda.bw && rm $output_dir/$sample/$sample"."_control_lambda.bdg\n";
    print SLURM "bedGraphToBigWig $output_dir/$sample/$sample"."_treat_pileup.bdg $chromsizes $output_dir/$sample/$sample"."_treat_pileup.bw && rm $output_dir/$sample/$sample"."_treat_pileup.bdg\n";
  
  close SLURM;
    system("sbatch $tmp_file");
  # unlink $tmp_file;
  
}

close(INPUT_FILE);
