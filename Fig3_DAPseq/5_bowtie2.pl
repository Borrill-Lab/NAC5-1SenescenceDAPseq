#!/usr/bin/perl -w

# 
#
# Aim of script is to run bowtie2 on DAP-seq data

#### paths and references:
my $path = '/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/bowtie2_index';
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/trimmed";
my $ref = "$path/iwgsc_refseq_1.1";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/mapped";

### lists of samples (text file containing directory/subdirectory with .fastq to map e.g. each line should look like: ERP004505/ERR392073/ in these subdirectories are the fastq.gz - text file must be in $output_dir):
#####
# Your list of samples should look like this:
# Sample_name read1 read2 

# Make sure the columns are tab separated


my $input_list_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/trimmed";
my $input_for_bowtie = "sample_list_trimmed_repeats.txt";


#######OUTPUT DIRECTORY, TMP DIRECTORY, SLURM_OUTPUT DIRECTORY MUST BE CREATED BEFORE RUNNING THE SCRIPT


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$input_list_dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_bowtie") || die "couldn't open the input file $input_for_bowtie!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
print "\nmy line was: $line\n";
  
#print "\nmy array: @array\n";
#print "\narray element 1: @array[0]\n";
  
my $sample = $array[0];
my $f1 = $array[1];
my $f2 = $array[2];
  
  
chdir("$read_path_triticum") or die "couldn't move to specific read directory $read_path_triticum";

#check if sample has already been processed
my $output_file = "$output_dir/$sample/$sample".".dupmark.sorted.bam.csi";
if (-e $output_file) {
        print "\nQuitting: output file $output_file already exists\n";
        next
}
  
  
my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel bowtie2 tasks
#
#SBATCH -p nbi-medium
#SBATCH -t 2-00:00
#SBATCH -c 4
#SBATCH --constraint=centos7
#SBATCH --mem=120000
#SBATCH -J bowtie2
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/mapped/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/mapped/slurm_output/%x.%N.%j.err
SLURM
  
 my $tmp_file = "$output_dir/tmp/bowtie2.$sample";
  
  
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";
  
  print SLURM "set -e\n";

  print SLURM "mkdir -p $output_dir/$sample\n";

  #bowtie2 2.4.0
  print SLURM "source package 29a74b59-88fc-4453-a30b-1310b34910b9\n";
  
    ### this part all needs editing!!!! ###
    print SLURM "bowtie2 --phred33 -q --fr --seed 55 -p 4 --rg-id $sample --rg 'SM:$sample' --rg 'LB:$sample' -x $ref -1 $f1 -2 $f2 -S $output_dir/$sample/$sample".".sam\n";

    #samtools 1.10
    print SLURM "source package aeee87c4-1923-4732-aca2-f2aff23580cc\n";

    print SLURM "samtools sort -@ 4 -T $output_dir/$sample/$sample".".temp_sort.bam -O bam -o $output_dir/$sample/$sample".".sorted.bam $output_dir/$sample/$sample".".sam\n";
    print SLURM ">&2 echo first sort\n";
    print SLURM "samtools index -c $output_dir/$sample/$sample".".sorted.bam\n";
    print SLURM ">&2 echo first index\n";
    print SLURM "samtools view -@ 4 -F 772 -q 30 -o $output_dir/$sample/$sample".".q.bam $output_dir/$sample/$sample".".sorted.bam\n";
    print SLURM ">&2 echo quality filtered\n";
    print SLURM "samtools sort -@ 4 -n -T $output_dir/$sample/$sample".".temp_namesort.bam -O bam -o $output_dir/$sample/$sample".".namesorted.bam $output_dir/$sample/$sample".".q.bam\n";
    print SLURM "samtools fixmate -@ 4 -m $output_dir/$sample/$sample".".namesorted.bam $output_dir/$sample/$sample".".fixmate.bam\n";
    print SLURM "samtools sort -@ 4 -T $output_dir/$sample/$sample".".qtemp_sort.bam -O bam -o $output_dir/$sample/$sample".".qsorted.bam $output_dir/$sample/$sample".".fixmate.bam\n";
    print SLURM ">&2 echo second sort\n";

    print SLURM "samtools markdup -@ 4 -d 2500 -f $output_dir/$sample/$sample".".markdup_metrics.txt $output_dir/$sample/$sample".".qsorted.bam -T $output_dir/$sample/$sample".".temp_markdup.bam $output_dir/$sample/$sample".".dupmark.bam\n";
    print SLURM ">&2 echo duplicates have been marked\n";
    print SLURM "samtools sort -@ 4 -T $output_dir/$sample/$sample".".duptemp_namesort.bam -O bam -o $output_dir/$sample/$sample".".dupmark.sorted.bam $output_dir/$sample/$sample".".dupmark.bam\n";
    print SLURM ">&2 echo third sort\n";
    print SLURM "samtools index -c $output_dir/$sample/$sample".".dupmark.sorted.bam\n";
    print SLURM ">&2 echo third index\n";

  close SLURM;
    system("sbatch $tmp_file");
  # unlink $tmp_file;
  
}

close(INPUT_FILE);
