#!/usr/bin/perl -w

# 
#
# Aim of script is to run macs2 on DAP-seq data

#### paths and references:
my $path = '/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/bowtie2_index';
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/mapped";
my $ref = "$path/iwgsc_refseq_1.1";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/peaks";

### lists of samples (text file containing directory/subdirectory with .fastq to map e.g. each line should look like: ERP004505/ERR392073/ in these subdirectories are the fastq.gz - text file must be in $output_dir):
#####
# Your list of samples should look like this:
# Sample_name	readpath

# Make sure the columns are tab separated


my $input_list_dir = "/jic/scratch/groups/Philippa-Borrill/Catherine/DAPseqPeaks/mapped";
my $input_for_macs2 = "input_list_macs2_52to54.txt";


#######OUTPUT DIRECTORY, TMP DIRECTORY, SLURM_OUTPUT DIRECTORY MUST BE CREATED BEFORE RUNNING THE SCRIPT


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$input_list_dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_macs2") || die "couldn't open the input file $input_for_macs2!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
print "\nmy line was: $line\n";
  
print "\nmy array: @array\n";
print "\narray element 1: @array[0]\n";
  
my $sample = $array[0];
my $f1 = $array[1];
#my $f2 = $array[2];
  
  
chdir("$read_path_triticum") or die "couldn't move to specific read directory $read_path_triticum";
  
  
my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel macs2 tasks
#
#SBATCH -p nbi-medium
#SBATCH -t 2-00:00
#SBATCH -c 1
#SBATCH --mem=120000
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

  #MACS 2.2.7.1
  print SLURM "source package 0aea19f1-78a3-40ee-87d3-49ab33c04ff5\n";
  
    ### this part all needs editing!!!! ###
    print SLURM "macs2 callpeak -t $f1 -f BAMPE --keep-dup auto --nomodel --broad --nolambda -g 12365172330 -n $sample --outdir $output_dir/$sample"."\n";

  
  close SLURM;
    system("sbatch $tmp_file");
  # unlink $tmp_file;
  
}

close(INPUT_FILE);
