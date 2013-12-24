#!/usr/bin/home python

########################################################################
#                                                                      #
#                           Alignment                                  #
#                                                                      #
#     This script aligns clean reads using novoalign.                  #
#     This script also allows for coverage files if desired.           #
#                                                                      #
#     Exteral dependencies: novoalign, samtools                        #
#                                                                      #
#     Script written by: Joshua Penalba (joshua.penalba@anu.edu.au)    #
#     Written on: 13 Aug 2013           Last Modified: 23 Dec 2013     #
#                                                                      #
########################################################################

import argparse
import sys
import os

aligner = argparse.ArgumentParser(description='Calls novoalign to align libraries')

aligner.add_argument('-f', dest='reads', help='directory of clean reads.')
aligner.add_argument('-r', dest='ref', help='directory of references.')
aligner.add_argument('-o', dest='out', type=str, help='output directory')
aligner.add_argument('-S', dest='maxscore', type=str, help='maximum score for alignment [90]', default = '90')
aligner.add_argument('-Z', dest='insize', type=str, help='library insert size [270]', default = '270')
aligner.add_argument('-D', dest='instd', type=str, help='library insert standard [70]', default = '70')
aligner.add_argument('-n', dest='maxlen', type=str, help='discard contigs longer than this (twice your read length) [200]', default='200')
aligner.add_argument('-b', dest='infofile', type=str, help='path to library info file (if batch processing is desired)')
aligner.add_argument('-cov', dest='cov', action='store_true', help='use flag if coverage files are desired')

if len(sys.argv)==1:
    aligner.print_help()
    sys.exit(1)
args = aligner.parse_args()
if args.infofile != None: args.homedir = '/'.join(args.infofile.split('/')[:-1])

if args.outdir.endswith('/'): pass
else: args.outdir = args.outdir+'/'

insertSize = args.insize
insertStd = args.instd
maxScore = args.maxscore
maxLen = args.maxlen

proj = set()
if args.infofile != None: 
    libinfo = open(args.infofile, 'r')
    for lines in libinfo:
        info = lines.strip().split()
        proj.add(args.homedir+info[2]+'/CleanedReads/')
    libinfo.close()
else: proj.add(args.libdir)

proj = list(proj)

for project in proj:
    if args.outdir != None: outdir = args.outdir
    else: outdir = '/'.join(project.split('/')[:-1])+'/Alignments/'
    try: os.mkdir(outdir)
    except OSError: pass
    if args.cov: 
        covpath = outdir+'Coverages/')
        try: os.mkdir(covpath)
        except OSError: pass
    if args.refpath != None: refpath = args.refpath
    else: refpath = '/'.join(project.split('/')[:-1])+'/References/'
    liblist = []
    allfiles = os.listdir(project)
    for eachfile in allfiles:
        if '_1_final.txt' in eachfile: liblist.append(project+eachfile)
        if eachfile.endswith('.gz'): os.system('gunzip %s%s' % (project, eachfile))
        else: pass
    for library in liblist:
        libname = library.split('/')[-1].split('_')[0]
        out = outdir+libname+'.sorted'
        read1 = library.replace('.gz', '')
        read2 = read1.replace('_1_', '_2_')
        readu = read1.replace('_1_', '_u_')
        ref = refpath+libname+'.fa'
#Call to index references
        os.system('novoindex novoindexref.out %s' % ref)

#Call to align paired reads
        os.system('novoalign -d novoindexref.out -f %s %s -i PE %s, %s -t %s -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -rAll -F STDFQ -o SAM > %soutPairedSam1' % (read1, read2, insertSize, insertStd, maxScore, outdir))

#Call to align unpaired reads
        os.system('novoalign -d novoindexref.out -f %s -t %s -n %s -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -rAll -F STDFQ -o SAM > %soutSoloSam1' % (readu, maxScore, maxLen, outdir))

#Grep only the aligned reads out of the novoalign output, then use samtools to merge and index
#only printout aligned reads
        os.system("grep -v \'ZS:Z:NM\' %soutPairedSam1 > %starget_pair.sam" % (outdir, outdir))
        os.system("grep -v \'ZS:Z:NM\' %soutSoloSam1 > %starget_solo.sam" % (outdir, outdir))
    
#run samtools to bam, merge sort and index
        os.system("samtools view -bS %starget_pair.sam > %starget_pair.bam" % (outdir, outdir))
        os.system("samtools view -bS %starget_solo.sam > %starget_solo.bam" % (outdir, outdir))
        os.system("samtools merge -f %starget.bam %starget_solo.bam %starget_pair.bam" % (outdir, outdir, outdir))
        os.system("samtools sort %starget.bam %s" % (outdir,out))
        os.system("samtools index %s.bam" % out)
        if args.cov: os.system("samtools depth %s.bam > %s%s.cov" % (out, covpath, library))
        os.system("gzip "+read1)
        os.system("gzip "+read2)
        os.system("gzip "+readu)
        os.system("rm %starget_pair.bam %starget_solo.bam %snovoindexref.out %soutSoloSam1 %soutPairedSam1 %starget.bam %starget_pair.sam %starget_solo.sam" % (outdir,outdir,outdir,outdir,outdir,outdir,outdir,outdir))

