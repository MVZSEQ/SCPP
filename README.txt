Sequence Capture using PCR-generated Probes (SCPP)
=================================================

Pipeline for Illumina SCPP experiment data.
Joshua Penalba, Sonal Singhal, & Ke Bi
December 24 2013

Pipeline is a constant work in progress. Contact Josh for python scripts and Ke for perl script questions.

This pipeline is a combination of novel scripts and scripts from other pipelines http://www.github.com/MVZSEQ/transcriptome and http://www.github.com/singhal/exomeCapture.

****************************************************************
1-pre-cleanup.pl: Reformats raw sequencing reads from Illumina for 2-scrubReads.pl

Dependencies: fastqc (if evaluation is desired)

Input is the raw libraries with the following naming convention: 
LIB1_ATCAGC_L006_R1_001.fasta.gz

Output in a folder called "pre-clean" will be:
LIB1_R1.fq

****************************************************************
2-scrubReads.pl: Removes adapters, contamination, and low complexity reads and merges overlapping reads

Dependencies: COPE, cutadapt, FLASH, fastQC, Bowtie2

Input is reformated libraries:
LIB1_R1.fq

Output is the following .txt files:
LIB1_1_final.txt (left reads)
LIB1_2_final.txt (right reads)
LIB1_u_final.txt (merged or unpaired reads)
LIB1.contam.out
LIB1.duplicates.out
LIB1.lowComplexity.out

Note: If only single indexing (P7) leave out P5 column

****************************************************************
