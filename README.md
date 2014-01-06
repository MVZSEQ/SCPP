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
LIB1_R2.fq

Output is the following .txt files:
LIB1_1_final.txt (left reads)
LIB1_2_final.txt (right reads)
LIB1_u_final.txt (merged or unpaired reads)
LIB1.contam.out
LIB1.duplicates.out
LIB1.lowComplexity.out

Note: If only single indexing (P7) leave out P5 column

****************************************************************
3-denovoAssembly.py: Uses ABySS to assemble each library.

Dependencies: ABySS, samtools

Input are cleaned reads:
LIB1_1_final.txt
LIB1_2_final.txt
LIB1_u_final.txt

Output is ABySS output:
LIB1_k61_5-contigs.fa
(k61 being the kmer and 5 being the kmer coverage)

****************************************************************
4-poolAssembly.py: Pools all the kmers and kmer coverages into one file/

Dependencies: none

Input are ABySS assemblies:
LIB1_k61_5-contigs.fa

Output:
LIB1.fa

****************************************************************
5-finalAssembly.pl: clusters contigs and assembles larger contigs.

Dependencies: blat (v34), CAP3, cd-hit-est

Input is the pooled assemblies:
LIB1.fa

Output:
LIB1.fa.final (final output)
LIB1.fa.original
LIB1.fa.cl98a
LIB1.fa.cl98b
LIB1.fa.cl98b.clustered
LIB1.fa.cl99a
LIB1.fa.cl99b

****************************************************************
6-generateReferences.py: generates a reference for each library using its own contigs.

Dependencies: bowtie2, samtools

Input:
LIB1_1_final.txt
LIB1.fa.final
RefLoci.fa

Output:
LIB1.fa
LIB1_ref.report

****************************************************************
7-alignment.py: Reads are aligned to the references using novoalign

Dependencies: novoalign, samtools

Input are cleaned reads and reference:
LIB1_1_final.txt
LIB1_2_final.txt
LIB1_u_final.txt
LIB1.fa

Output:
LIB1.sorted.bam
LIB1.sorted.bam.bai
LIB1.cov (only if coverage is desired)

****************************************************************
8-phaseHap.pl: Uses GATK to come up with haplotypes

Dependencies: GATK, picard

Input are alignments:
LIB1.sorted.bam
LIB1.sorted.bam.bai
LIB1.fa (reference file)

Output: 
LIB1.final_hap1.txt
LIB1.final_hap2.txt

****************************************************************
9-hapTools.py: Has different functions to manipulate the haplotypes

Command: callFixer
Function: fixes the references so that each base is the majority of all the reads that map. This is how homozygous loci are called. The mitochondrial haplotype is the output of this script. The script needs to be called for nuclear loci as well so that the homozygous loci are not comprised of bases from reads with sequencing errors. 

Dependencies: samtools, assembly_evaluation.pl

Input are alignments and references:
LIB1.fa (reference file)
LIB1.sorted.bam
LIB1.sorted.bam.bai

Output:
LIB1_fixed.fa
LIB1_pos.txt


Command: consensus
Function: Creates a consensus sequences from the haplotypes. It is always preferable to use the haplotype data from GATK rather than consensus sequences. This is to be used by the researcher's discretion. The output is only the heterozygous loci. The output should be combined with fixed references using the 'combine' command to get all the loci.

Dependecies: mafft

Input:
LIB1.final_hap1.txt
LIB1.final_hap2.txt

Output: 
LIB1.fa (heterozygous consensus)
LIB1.out (output information with indel and SNP info)


Command: combine
Function: Combines the heterozygous output of the 'consensus' function with the fixed references to get all the loci in one file.

Dependencies: none

Input:
LIB1.fa (heterozygous only)

Output:
LIB1.fa (combined)


Command: rename
Function: Renames the libraries if the name is not what the researcher desires.

Input:
LIB1.fa (combined or otherwise)
ExtInfo.txt

Output:
RENAME1.fa


Command: perLoc
Function: Turns each library consensus file to files organized per locus rather than per library.

Input:
LIB1.fa (combined or RENAME1.fa)

Output:
LOC1.fa
LOC2.fa
LOC3.fa

