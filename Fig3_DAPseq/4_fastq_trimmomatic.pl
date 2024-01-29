#!/usr/bin/perl -w

# 
#
# Aim of script is to run trimmomatic on DAP-seq data
# and generate FASTQC reports

#### paths and references:
my $path = '/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/bowtie2_index';
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/raw";
my $ref = "$path/iwgsc_refseq_1.1_parts";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/trimmed";
my $fastqc_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/data/FASTQC";

### lists of samples (text file containing directory/subdirectory with .fastq to map e.g. each line should look like: ERP004505/ERR392073/ in these subdirectories are the fastq.gz - text file must be in $input_list_dir):
#####
# Your list of samples should look like this:
# Sample_name	read1	read2 

# Make sure the columns are tab separated


my $input_list_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/raw";
my $input_for_trimmomatic = "sample_list_53and54.txt";


#######OUTPUT DIRECTORY, TMP DIRECTORY, SLURM_OUTPUT DIRECTORY MUST BE CREATED BEFORE RUNNING THE SCRIPT


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$input_list_dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_trimmomatic") || die "couldn't open the input file $input_for_trimmomatic !";
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
my $output_file = "$output_dir/$sample/$sample"."_P1.fastq.gz";
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
#SBATCH -t 0-03:00
#SBATCH -c 4
#SBATCH --mem=30000
#SBATCH -J trimmomatic
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/trimmed/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/fastq/trimmed/slurm_output/%x.%N.%j.err
SLURM
  
 my $tmp_file = "$output_dir/tmp/trimmomatic.$sample";
  
  
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";
  
  print SLURM "set -e\n";

    ### this part all needs editing!!!! ###

    #FastQC 0.11.9
    print SLURM "source package a138abc4-ed60-4774-8db8-80b4770b1710\n";
    #Trimmomatic 0.39
    print SLURM "source package 50fcf79b-73a3-4f94-9553-5ed917823423\n";

    print SLURM "mkdir -p $output_dir/$sample\n";
    print SLURM "mkdir -p $fastqc_dir/trimmed/$sample\n";

    print SLURM "fastqc -t 2 -o $fastqc_dir/raw $f1 $f2\n";

    print SLURM "trimmomatic PE -threads 4 -phred33 $f1 $f2"." \\\n";
    print SLURM "-baseout $output_dir/$sample/$sample".".fastq.gz \\\n";
    print SLURM "ILLUMINACLIP:/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33\n"; 

    print SLURM "fastqc -t 4 -o $fastqc_dir/trimmed/$sample $output_dir/$sample/$sample"."_1P.fastq.gz $output_dir/$sample/$sample"."_2P.fastq.gz $output_dir/$sample/$sample"."_1U.fastq.gz $output_dir/$sample/$sample"."_2U.fastq.gz\n";
  
  close SLURM;
    system("sbatch $tmp_file");
  # unlink $tmp_file;
  
}

close(INPUT_FILE);
