#!/usr/bin/home python

#############################################################################
#                                                                           #
#                        Generating Reference                               #
#                                                                           #
#  This script is to generate references for mapping PCR amplicon bait      #
#  capture data. It is only necesarry to do this if there is no available   #
#  sequence data of the target loci from closely related individuals.       #
#                                                                           #
#  External dependencies: bowtie2, samtools, BLAST                          #
#                                                                           #
#  Script written by: Joshua Penalba (joshua.penalba@anu.edu.au)            #   
#  Written on: 27 Jul 2013                Last modified: 23 Dec 2013        #                                                           
#                                                                           # 
#############################################################################

import argparse
import os
import sys

ref = argparse.ArgumentParser(description='Recovers final assembled contigs that corresponds to each target for a SCPP project')

ref.add_argument('-f', dest='reads', type=str, help='directory of cleaned reads.')
ref.add_argument('-r', dest='ref', type=str, help='fasta file of target sequences', required = True)
ref.add_argument('-a', dest='assembly', type=str, help='directory of final assemblies (ends in .fasta.final)')
ref.add_argument('-o', dest='outdir', type=str, help='output directory')
ref.add_argument('-d', dest='readhead', type=str, help='beginning of the read header (ex. @M00384, @HWI, @HS) [@HWI]', default = '@HWI')
ref.add_argument('-e', dest='evalue', type=str, help='e-value for BLAST search [1e-10]', default = '1e-10')
ref.add_argument('-b', dest='infofile', type=str, help='library info file. Only use if batch processing is desired. Assumes standard SCPP directory nomenclature')
ref.add_argument('-nu', dest='nuc', action="store_true", help='loci are all only nuclear loci (only if batch processing)')
ref.add_argument('-mt', dest = 'mito', action="store_true", help='loci are all mitochondrial (only if batch processing)')

if len(sys.argv) == 1:
    ref.print_help()
    sys.exit(1)
args = ref.parse_args()

if args.outdir.endswith('/'): pass
else: args.outdir = args.outdir+'/'
if args.infofile != None:
    if args.dir != None or args.outdir != None or args.assembly != None: 
        print 'If batch processing by info file is selected, flags -f, -o, and -a should not be activated.'
        sys.exit()
elif args.infofile == None:
    if args.dir == None or args.outdir == None or args.assembly == None:
        print 'If single directory proceessing is desired, flags -f, -o, and -a should ALL be assigned a path.'
        sys.exit()

if args.infofile: args.homedir = '/'.join(args.infofile.split('/')[:-1])
if args.nuc: appname = '_nu'
elif args.mito: appname = '_mt'
else: appname = ''

######################
# ESTABLISHING PATHS #
######################

liblist = set()
libs = {}

if args.dir != None:
    allfiles = os.listdir(args.dir)
    for each in allfiles:
        if '_final.txt' in each:
            library = each.split('_')[0]
            liblist.add(library)
elif args.infofile != None:
    libinfo = open(args.infofile, 'r')
    for lines in libinfo:
        info = lines.strip().split()
        liblist.add(info[0])
        libs[info[0]] = info[2] 

evalue = args.evalue
readhead = args.readhead
ToTargets = args.ref
liblist = list(liblist)


for lib in liblist:
    if args.infofile != None and args.dir == None:
        proj = libs[lib]
        reads1 = '%s%s/CleanedReads/%s_1_final.txt' % (args.homedir, proj, lib)
        assembly = '%s/%s/FinalAssemblies/%s.fa.final' % (args.homedir, proj, lib)
        OutDir = '%s/%s/References%s' % (args.homedir, proj, appname)
    elif args.dir != None and args.infofile == None:
        reads1 = '%s%s_1_final.txt' % (args.dir, lib)
        OutDir = args.outdir
        assembly = '%s%s.fa.final' % (args.assembly, lib)
    print 'PROCESSING %s' % lib

    try: os.mkdir(OutDir)
    except: OSError
    pass

###################
# Initial Mapping #
###################

    print 'RUNNING INITIAL MAPPING...'

    TargPath = '/'.join(ToTargets.split('/')[:-1])
    if ToTargets+'.1.bt2' in os.listdir(TargPath): pass
    else: os.system("bowtie2-build -q %s %s" % (ToTargets, ToTargets))
    os.system('samtools faidx '+ToTargets)

    reads2 = reads1.replace('_1','_2')
    readsu = reads1.replace('_1','_u')
    os.system('gunzip '+reads1+'.gz')
    os.system('gunzip '+reads2+'.gz')
    os.system('gunzip '+readsu+'.gz')
    #lib = '_'.join(reads1.split('/')[-1].split('_')[0:2])

    #os.system("bowtie2 -x %s -1 %s -2 %s -S bowtie1.sam -5 5 -3 5 --very-sensitive-local -k 10 -X 300 -p 4" % (ToTargets, reads1, reads2))
    #os.system("bowtie2 -x %s %s -S bowtie2.sam -5 5 -3 5 --very-sensitive-local -k 10 -p 4" % (ToTargets, readsu))
    os.system("bowtie2 -x %s -1 %s -2 %s -S %s/bowtie1.sam -5 5 -3 5 --very-sensitive-local -k 10 -X 300 -p 4" % (ToTargets, reads1, reads2, OutDir))
    os.system("bowtie2 -x %s %s -S %s/bowtie2.sam -5 5 -3 5 --very-sensitive-local -k 10 -p 4" % (ToTargets, readsu, OutDir))
    print 'MAPPING COMPLETE. RUNNING SAMTOOLS FILE CONVERSIONS...'
    os.system("samtools view -bS %s/bowtie1.sam > %s/bowtie1.bam" % (OutDir, OutDir))
    os.system("samtools view -bS %s/bowtie2.sam > %s/bowtie2.bam" % (OutDir, OutDir))
    os.system("samtools merge %s/bowtie.bam %s/bowtie1.bam %s/bowtie2.bam" % (OutDir, OutDir, OutDir))
    os.system("samtools sort %s/bowtie.bam %s/bowtie.sorted" % (OutDir, OutDir))
    os.system("rm %s/bowtie1.sam %s/bowtie2.sam %s/bowtie.bam %s/bowtie1.bam %s/bowtie2.bam" % (OutDir, OutDir, OutDir, OutDir, OutDir))
    os.system("samtools index %s/bowtie.sorted.bam" % OutDir)
    os.system("samtools view %s/bowtie.sorted.bam > %s/bowtie.sam" % (OutDir, OutDir))
    os.system("rm %s.* %s/bowtie.sorted*" % (ToTargets, OutDir))

######################
# Parse Mapped Reads #
######################

    print 'SAMTOOLS FILE CONVERSIONS COMPLETE. PARSING MAPPED READS...'

    samfile = open(OutDir+'/bowtie.sam', 'r')
    reads_1 = open(reads1, 'r')
    readsout = open(OutDir+'/readsout.fa', 'w')
    refloci = open(ToTargets, 'r')
    readlib = {}

#Creates a library with all read that mapped to each locus
    for lines in samfile:
        info = lines.strip().split()
        if info[2] not in readlib:
            readlib[info[2]] = ['@'+info[0]]
        elif info[2] in readlib:
            readlib[info[2]].append('@'+info[0])
        else: pass
    samfile.close()

# Create a list of all the loci
    allloci = []
    for lines in refloci:
        info = lines.strip().split()
        if len(info) == 1 and info[0][0] == '>': allloci.append(info[0][1:])
        else: pass
    allloci.sort()
    refloci.close()
    refloci = open(ToTargets , 'r')
    reflib = {}
    for lines in refloci:
        info = lines.strip().split()
        if len(info) == 0: continue
        elif info[0].startswith('>'): locus = info[0][1:]
        elif locus not in reflib: reflib[locus] = [info[0]]
        elif locus in reflib: reflib[locus].append(info[0])
    for contig in reflib:
        if len(reflib[contig]) == 1: reflib[contig] = reflib[contig][0]
        elif len(reflib[contig]) > 1: reflib[contig] = ''.join(reflib[contig])
    refloci.close()

#Create a list of all the loci that are not found
    notfound = []
    for keys in readlib:
        if keys in allloci: continue
        else: notfound.append(keys)
    for things in notfound:
        del readlib[things]

# Create a library for output report
    report = {}
    for locus in allloci: 
        if locus in readlib: report[locus] = ['YES']
        else: report[locus] = ['NO']

#Create a library with contigs as keys instead
    readlib2 = {}
    for keys in readlib:
        for reads in readlib[keys]:
            readlib2[reads] = keys

#Output for BLAST
    recording = 'OFF'
    for lines in reads_1:
        info = lines.strip().split()
        if info[0].startswith(readhead) and info[0].split('/')[0] in readlib2:
            readsout.write('>'+readlib2[info[0].split('/')[0]]+'\n')
            recording = 'ON'
        elif info[0].startswith('+'): 
            recording = 'OFF'
            continue
        elif recording == 'ON':
            readsout.write(info[0]+'\n')
        else: pass
    readsout.close()
    reads_1.close()

############################################
# BLAST mapped reads onto final assemblies #
############################################

    print ('PARSING READS COMPLETE. BLAST TO FINAL ASSEMBLY...')

    os.system("formatdb -i %s -p F" % (assembly))
    os.system("blastall -p blastn -d %s -i %s/readsout.fa -a 4 -e %s -m 8 -o %s/blast.out -b 10" % (assembly, OutDir, evalue, OutDir))
    os.system("rm %s.*" % (assembly))

#########################################
# Generating final reference and report #
#########################################

    print ('BLAST COMPLETE. GENERATING FINAL REFERENCE...')

    blastout = open(OutDir+'/blast.out','r')
    contigs = open(assembly, 'r')
    rep = open(OutDir+'/'+lib+'_ref.report', 'w')

#Library of BLAST hits
    blastlib = {}
    for lines in blastout:
        info = lines.strip().split()
        if info[0] not in blastlib:
            blastlib[info[0]] = [info[1]]
        elif info[0] in blastlib:
            blastlib[info[0]].append(info[1])
    blastout.close()

#Set of contighits
    contighits = set()
    for locus in blastlib:
        for read in blastlib[locus]:
            contighits.add(read)

#Library of contigs
    contiglib = {}
    for lines in contigs:
        info = lines.strip().split()
        if info[0].startswith('>') and info[0][1:] in contighits:
            contiglib[info[0][1:]] = []
            contigname = info[0][1:]
            recording = 'ON'
        elif info[0].startswith('>') and info[0][1:] not in contighits:
            recording = 'OFF'
        elif recording == 'ON': contiglib[contigname].append(info[0])
        else: pass
    for contig in contiglib:
        if len(contiglib[contig]) == 1: contiglib[contig] = contiglib[contig][0]
        elif len(contiglib[contig]) > 1: contiglib[contig] = ''.join(contiglib[contig])
#Library of contig lengths
    contiglenlib = {}
    for contig in contiglib:
        length = len(contiglib[contig])
        contiglenlib[contig] = length

#Library of contig hits per locus
    mosthits = {}
    for locus in blastlib:
        mosthits[locus] = {}
        for read in blastlib[locus]:
            if read not in mosthits[locus]: mosthits[locus][read] = 1
            elif read in mosthits[locus]: mosthits[locus][read] += 1
#Library of all top hits
    tophits = {}
    for locus in mosthits:
        numhit = 0
        for read in mosthits[locus]:
            hits = mosthits[locus][read]
            if hits > numhit: numhit = hits
            else: continue
        choices = []
        for read1 in mosthits[locus]:
            hits1 = mosthits[locus][read1]
            if hits1 == numhit: choices.append(read1)
            else: continue
        longest = 0
        for read2 in choices:
            if contiglenlib[read2] > longest:
                longest = contiglenlib[read2]
                tophit = read2
            else: continue
        tophits[locus] = tophit
        report[locus].append(tophit)

#Generating final file
    refout = open(OutDir+'/'+lib+'.fa','w')
    for locus in tophits:
        refout.write('>'+locus+'\n')
        refout.write(contiglib[tophits[locus]]+'\n')


#Generating Report
    rep.write('Locus\tMapped Initially?\tFinal Reference\n')
    for locus in allloci:
        if len(report[locus]) == 2: rep.write(locus+'\t'+'\t'.join(report[locus])+'\n')
        else: 
            rep.write(locus+'\t'+'\t'.join(report[locus])+'\tNot in final\n')
        #refout.write('>%s\n' % locus)
        #refout.write(reflib[locus]+'\n')
    rep.close()
    refout.close()
    os.system('gzip %s' % reads1)
    os.system('gzip %s' % reads2)
    os.system('gzip %s' % readsu)
    os.system("rm %s/readsout.fa %s/blast.out" % (OutDir, OutDir))
    os.system("rm %s/bowtie.sam" % (OutDir))
    print('FINAL REFERENCE AND REPORT COMPLETED.')
