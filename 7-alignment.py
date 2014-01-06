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
#     Written on: 13 Aug 2013           Last Modified: 3 Jan 2014      #
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
aligner.add_argument('-cov', dest='cov', action='store_true', help='use flag if coverage files are desired')

if len(sys.argv)==1:
    aligner.print_help()
    sys.exit(1)
args = aligner.parse_args()

if args.out.endswith('/'): pass
else: args.out = args.out+'/'

insertSize = args.insize
insertStd = args.instd
maxScore = args.maxscore
maxLen = args.maxlen

proj = set()
proj.add(args.reads)
proj = list(proj)

for project in proj:
    outdir = args.out
    try: os.mkdir(outdir)
    except OSError: pass
    try: os.mkdir(outdir+'tmp')
    except OSError: pass
    if args.cov: 
        covpath = outdir+'Coverages/'
        try: os.mkdir(covpath)
        except OSError: pass
    refpath = args.ref
    liblist = []
    allfiles = os.listdir(project)
    for eachfile in allfiles:
        if '_1_final.txt' in eachfile: liblist.append(project+eachfile)
        if eachfile.endswith('.gz'): os.system('gunzip %s%s' % (project, eachfile))
        else: pass
    for library in liblist:
        libname = library.split('/')[-1].split('_')[0]
        out = outdir+libname+'.sorted'
        read1 = outdir+'tmp/'+libname+'_1_final.txt'
        read2 = read1.replace('_1_', '_2_')
        readu = read1.replace('_1_', '_u_')
        os.system('cp %s* %stmp/' % (library.split('1_final')[0], outdir))
        ref = refpath+libname+'.fa'
#Call to index references
        os.system('novoindex %s%s.out %s' % (outdir, libname, ref))

#Call to align paired reads
        os.system('novoalign -d %s%s.out -f %s %s -i PE %s, %s -t %s -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -rAll -F STDFQ -o SAM > %soutPairedSam1' % (outdir, libname, read1, read2, insertSize, insertStd, maxScore, outdir))

#Call to align unpaired reads
        os.system('novoalign -d %s%s.out -f %s -t %s -n %s -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -rAll -F STDFQ -o SAM > %soutSoloSam1' % (outdir, libname, readu, maxScore, maxLen, outdir))

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
        if args.cov: os.system("samtools depth %s.bam > %s%s.cov" % (out, covpath, libname))
        else: pass
        os.system("gzip "+read1)
        os.system("gzip "+read2)
        os.system("gzip "+readu)
        os.system("rm %starget_pair.bam %starget_solo.bam %s%s.out %soutSoloSam1 %soutPairedSam1 %starget.bam %starget_pair.sam %starget_solo.sam %stmp/*" % (outdir,outdir,outdir,libname,outdir,outdir,outdir,outdir,outdir, outdir))
    os.system('rmdir %stmp' % outdir)
