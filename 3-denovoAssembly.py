#!/usr/bin/home/python

##############################################################################
#                                                                            #
#                    De novo Assembly of the Samples                         #
#                                                                            #
#   This script does an ABySS assembly of each species for a SCPP            #
#   type project.                                                            #
#                                                                            #
#   External dependencies: ABySS, samtools                                   #
#                                                                            #
#   Script written by: Joshua Penalba (joshua.penalba@anu.edu.au)            #
#   Written on: 13 Aug 13              Last Modification: 23 Dec 13          #
#                                                                            #
##############################################################################


import argparse
import sys
import os

assem = argparse.ArgumentParser(description='Performs a multi-kmer de novo assembly of the samples using ABySS')

assem.add_argument('-f', dest='reads', type=str, help='directory of cleaned reads.')
assem.add_argument('-o', dest='out', type=str, help='output directory.')
assem.add_argument('-k', dest='kmers', default = ['21', '41', '61', '81'], nargs='*', help='list of kmers. [21,41,61,81]')
assem.add_argument('-c', dest='kcov', default = ['5','10','20'], nargs='*', help='list of kmer coverages. [5,10,20]')
assem.add_argument('-b', dest='infofile', type=str, help='path to library info file. Only use if batch processing is desired. Assumes standard SCPP data organization.')
assem.add_argument('-p', dest='numproc', type=str, help='number of processors. [8]', default=str(8))

if len(sys.argv)==1:
    assem.print_help()
    sys.exit(1)
args = assem.parse_args()

if args.infofile: args.homedir = '/'.join(args.infofile.split('/')[:-1])
if not args.out.endswith('/'): args.out = args.out + '/'
if not args.reads.endswith('/'): args.reads = args.reads + '/'

##################
# ABYSS ASSEMBLY #
##################

#CREATES LIST OF INFO
proj = set()
if args.infofile != None:
    libinfo = open(args.infofile, 'r')
    for lines in libinfo:
        info=lines.strip().split()
        proj.add(args.homedir+info[2]+'/CleanedReads/')
    libinfo.close()
else: proj.add(args.reads)

proj = list(proj)
workdir = os.getcwd()

#CREATES FOLDERS AND RUNS ASSEMBLY
for project in proj:
    if args.out != None: outdir = args.out
    else: outdir = '/'.join(project.split('/')[:-1])+'/Assemblies/'
    try: os.mkdir(outdir)
    except OSError: pass
    liblist = []
    allfiles = os.listdir(project)
    for eachfile in allfiles:
        if '_1_final.txt' in eachfile: liblist.append(project+eachfile)
        if eachfile.endswith('.gz'): os.system('gunzip %s%s' % (project, eachfile))    
        else: pass
    for library in liblist:
        libname = library.split('/')[-1].split('_')[0]
        reads1 = library.replace('.gz', '')
        reads2 = reads1.replace('_1_', '_2_')
        readsu = reads1.replace('_1_', '_u_')
        for kmer in args.kmers:
            for kcov in args.kcov:
                os.system("abyss-pe k=%s np=%s n=5 s=200 in='%s %s' se='%s' name=%s%s_k%s_%s c=%s e=%s E=0" % (kmer, args.numproc, reads1, reads2, readsu, outdir, libname, kmer,kcov, kcov, kcov))
        os.system('gzip %s' % reads1)
        os.system('gzip %s' % reads2)
        os.system('gzip %s' % readsu)
        os.chdir(workdir)
    
os.system('rm *.log')
os.system('rm target*')
os.system('rm coverage.hist')
os.system('rm *Sam*')
