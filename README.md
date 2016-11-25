Travis CI Status: [![Build Status](https://travis-ci.org/serosko/BAMQC.svg?branch=master)](https://travis-ci.org/serosko/BAMQC)

BAMQC - Simple quality-control for (single-sample) BAM-files.
=============================================================

SYNOPSIS  

    BAMQC BAM_FILE [OPTIONS]

DESCRIPTION  

    -h, --help  
          Display the help message.
    --version  
          Display version information.

  I/O Options:  

    -r, --reference IN  
          Path to reference genome. Required for C>A/G>T-Artifact-check. 
	  Valid filetypes are: fasta, fa, fastq, fq, fasta.gz, fa.gz, fastq.gz, fq.gz, fasta.bz2, fa.bz2, fastq.bz2, and fq.bz2. 
 
    -oi, --output-file-inserts OUT  
          Path to output file for the insert-size distribution.

    -oc, --output-file-conversions OUT  
          Path to output file for the C>A/G>T-Artifact-check.

  General Options:  

    -mmq, --min-mapq INT  
          Minimum mapping quality. In range [0..244]. Default: 25.

    -v0, --no-verbosity  
          Disable parameter feedback.

  Insert-size-distribution Options:  

    -i, --insert-size-distribution  
          Counts the insert size of each valid read-pair. Output to standard output if -oi with path is not specified.

    -m, --max-insert INT  
          Maximum insert size. Sizes above will be ignored. In range [100..inf]. Default: 1000.

  C>A/G>T-Artifact Options:  

    -c, --cagt-artifact  
          Perform check for C>A/G>T artifacts induced during sample preparation (Costello et al. (2013)).Requires
          reference genome. Output to standard output if -oc with path is not specified.

EXAMPLES  

    BAMQC file.bam -i  
          Calculate insert-size distribution of mapped reads in "file.bam" and print output to standard output.

    BAMQC file.bam -r reference.fa -i -c  
          Check if reads in "file.bam" show an increased ammount of strand specific C>A/G>T conversions. and print
          output to standard output.

    BAMQC file.bam -r reference.fa -oi inserts.txt -oc conversions.txt -mmq 30 -m 500  
          Perform both of the above checks in one run, but only consider alignments with mapping-quality of at least
          phred 30 (both checks) and a maximum considered insert-size of 500 (insert-size calculation). Output is
          saved in "inserts.txt" and "conversions.txt."

VERSION  
    Last update: Nov 15 2016  
    BAMQC version: 1.0.0  
    SeqAn version: 2.2.0  
    Author: Sebastian Roskosch <Sebastian.Roskosch[at]bihealth.de>
