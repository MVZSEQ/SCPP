#!/usr/bin/home python

# This script is to deal with haplotypes

# SETTING UP ARGUMENT PARSING
###############################

import argparse
import sys
import os

hapTools = argparse.ArgumentParser(description = 'Various tools to deal with SCPP haplotype data')
subparsers = hapTools.add_subparsers(dest='command')

hapTools_cf = subparsers.add_parser('callFixer', help='Calls the assembly fixer script. [dependencies = samtools]')
hapTools_cf.add_argument('-a', dest='bam', type=str, help='directory of bam files')
hapTools_cf.add_argument('-r', dest='ref', type=str, help='directory of reference files')
hapTools_cf.add_argument('-o', dest='out', type=str, help='output directory')
hapTools_cf.add_argument('-s', dest='script', type=str, help='path to assembly_evaluation.pl script')

hapTools_con = subparsers.add_parser('consensus', help='Turns haplotypes into consensus sequences. [dependencies = mafft]')
hapTools_con.add_argument('-a', dest='hap', type=str, help='directory of haplotypes.')
hapTools_con.add_argument('-o', dest='out', type=str, help='output directory.')

hapTools_co = subparsers.add_parser('combine', help='Combines the homozygotes from the fixed references to the heterozygote output of this data')
hapTools_co.add_argument('-a', dest='con', type=str, help='directory of consensus sequences')
hapTools_co.add_argument('-f', dest='ref', type=str, help='directory of fixed references')
hapTools_co.add_argument('-o', dest='out', type=str, help='output directory')

hapTools_rn = subparsers.add_parser('rename', help='Renames the library names to extraction names')
hapTools_rn.add_argument('-a', dest='con', type=str, help='directory of consensus sequences')
hapTools_rn.add_argument('-s', dest='sampname', type=str, help='file with sample name.')
hapTools_rn.add_argument('-o', dest='out', type=str, help='output directory')

hapTools_pl = subparsers.add_parser('perLoc', help='Takes consensus sequences organized per library and converts them per locus with the name as the header')
hapTools_pl.add_argument('-f', dest='dir', type=str, help='path to directory of consensus sequences')
hapTools_pl.add_argument('-o', dest='outdir', type=str, help='output directory [default:same as input]')


if len(sys.argv)==1:
    hapTools.print_help()
    sys.exit(1)

args = hapTools.parse_args()

if len(sys.argv)==2:
    if args.command == 'consensus': hapTools_con.print_help()
    elif args.command == 'rename': hapTools_rn.print_help()
    elif args.command == 'perLoc': hapTools_pl.print_help()
    elif args.command == 'combine': hapTools_co.print_help()
    elif args.command =='callFixer': hapTools_cf.print_help()
    sys.exit(1)

#Command that uses IUPAC code to code for SNPs                                                                                                                                       
def SNP(x,y):
    z = x+y
    if z == 'AC' or z == 'CA': a = 'M'
    elif z == 'AG' or z == 'GA': a = 'R'
    elif z == 'AT' or z == 'TA': a = 'W'
    elif z == 'CG' or z == 'GC': a = 'S'
    elif z == 'CT' or z == 'TC': a = 'Y'
    elif z == 'GT' or z == 'TG': a = 'K'
    elif z.find('-') != -1:
        if z.find('A') != -1: a = 'M'
        elif z.find('G') != -1: a = 'M'
        elif z.find('T') != -1: a = 'M'
        elif z.find('C') != -1: a = 'M'
        elif z.find('N') != -1: a = ''
    elif z.find('N') != -1:
        if z.find('A') != -1: a = 'A'
        elif z.find('G') != -1: a = 'G'
        elif z.find('T') != -1: a = 'T'
        elif z.find('C') != -1: a = 'C'
    return  a

###############
#  CONSENSUS  #
###############

if args.command == 'consensus':
    if args.hap.endswith('/'): pass
    else: args.hap = args.hap+'/'
    if args.out.endswith('/'): pass
    else: args.out = args.out+'/'
    allfiles = os.listdir(args.hap)
    libs = []
    for each in allfiles:
        if each.endswith('_hap1.txt'): libs.append(each.split('.')[0])
        else: pass
    for library in libs:
        mastlib = {}
        uEvent = []
        for i in range(1,3):
            haplo = open(args.hap+'%s.final_hap%s.txt' % (library, i), 'r')
            haplib = {}
            vcflib = {}
            for lines in haplo:
                info = lines.strip().split()
                if info[0].startswith('>'):
                    if len(info[0].split('_')) > 1:
                        vcfhead = info[0].split('_')
                        locus = vcfhead[0][1:]
                        order = vcfhead[2]
                        lib = 'vcf'
                    else:
                        locus = info[0][1:]
                        lib = 'hap'
                else:
                    seq = info[0]
                    if lib == 'vcf' and locus not in vcflib:
                        vcflib[locus] = {order:seq.strip('N')}
                    elif lib == 'vcf' and locus in vcflib:
                        vcflib[locus][order] = seq.strip('N')
                    else: haplib[locus] = seq
            haplo.close()
#Merges the vcf files                                                                                                                                                               
            for locus in vcflib:
                ordering = vcflib[locus].keys()
                ordering.sort()
                tmpseq = []
                for each in ordering:
                    tmpseq.append(vcflib[locus][each])
                newseq = ''.join(tmpseq)
                haplib[locus] = newseq
            mastlib[i] = haplib

#Make temporary fasta files per locus
        try: os.mkdir(args.out)
        except OSError: pass
        try: os.mkdir(args.out+library)
        except OSError: pass
        os.chdir(args.out+library)

        for locus in mastlib[1]:
            locfile = open(locus+'.txt', 'w')
            seq1 = mastlib[1][locus]
            seq2 = mastlib[2][locus]
            if '<UNASSEMBLED_EVENT>' in seq1:
                uEvent.append(locus)
                continue
            elif '<UNASSEMBLED_EVENT>' in seq2:
                uEvent.append(locus)
                continue
            locfile.write('>%s\n%s\n>%s\n%s' % (locus, seq1, locus, seq2))
            locfile.close()
            os.system('mafft %s.txt > %s.fa' % (locus, locus))

#Finally makes the consensus files 
        conout = open(args.out+library+'.fa','w')
        conrep = open(args.out+library+'.out','w')
        for locus in mastlib[1]:
            mkey = 0
            try: mfile = open(locus+'.fa', 'r')
            except IOError: continue
            mseq = {}
            conseq = []
            for lines in mfile:
                info = lines.strip().split()
                if info[0].startswith('>'):
                    mkey += 1
                else:
                    if mkey in mseq: mseq[mkey].append(info[0])
                    else: mseq[mkey] = [info[0]]
            for i in range(1,3):
                try: mseq[i] = ''.join(mseq[i]).upper()
                except KeyError:
                    print locus+' is being read as empty'
                    sys.exit()
            for bp in enumerate(mseq[1]):
                if bp[1] == mseq[2][bp[0]]: conseq.append(bp[1])
                else:
                    if bp[1] or mseq[2][bp[0]] == '-':
                        conrep.write('%s\t%s\t%s%s\n' % (locus, bp[0], bp[1],mseq[2][bp[0]]))
                    conseq.append(SNP(bp[1],mseq[2][bp[0]]))
            conseq = ''.join(conseq)
            conout.write('>%s\n%s\n' % (locus, conseq))
        for locus in uEvent:
            conrep.write('%s has an unassembled event\n' % locus)
        conout.close()
        conrep.close()
        os.system('rm -r %s%s' % (args.out, library))

#############
#  COMBINE  #
#############

# Combines consensus sequences from this script (hets only) with the homozygotes from fixed references.

elif args.command == 'combine':
    if args.out.endswith('/'): pass
    else: args.out = args.out+'/'
    if args.con.endswith('/'): pass
    else: args.con = args.con+'/'
    if args.ref.endswith('/'): pass
    else: args.ref = args.ref+'/'
    conFiles = os.listdir(args.con)
    for each in conFiles:
        dataDir = {}
        if each.endswith('.fa'):
            lib = each.split('.')[0]
            consensus = open(args.con+each, 'r')
            for lines in consensus:
                info = lines.strip().split()
                if info[0].startswith('>'): 
                    locus = info[0][1:]
                    continue
                else: dataDir[locus] = info[0]
            consensus.close()
            fixed = open(args.ref+lib+'_fixed.fa', 'r')
            for lines in fixed:
                info = lines.strip().split()
                if info[0].startswith('>'):
                    locus = info[0][1:]
                    continue
                else:
                    if locus in dataDir: continue
                    else: dataDir[locus] = info[0]
            fixed.close()
        else: continue
        try: os.mkdir(args.out)
        except OSError: pass
        outfile = open(args.out+lib+'.fa', 'w')
        for locus in dataDir:
            outfile.write('>%s\n%s\n' % (locus, dataDir[locus]))
        outfile.close()

############
#  RENAME  #                                                                                                    
############
                        
elif args.command == 'rename':
    samplib = {}
    infofile = open(args.sampname, 'r')
    if args.out.endswith('/'): pass
    else: args.out = args.out+'/'
    try: os.mkdir(args.out)
    except OSError: pass
    for lines in infofile:
        info = lines.strip().split()
        samplib[info[0]] = info[1]
    infofile.close()
    allfiles = os.listdir(args.con)
    libfiles = []
    for each in allfiles:
        if each.endswith('.fa'): libfiles.append(each.split('.')[0])
        else: pass
    for library in libfiles:
        os.system('cp %s%s*fa %s%s.fa' % (args.con, library, args.out, samplib[library.split('_')[0]]))

############
#  PERLOC  #
############

elif args.command == 'perLoc':
    if args.outdir.endswith('/'): pass
    else: args.outdir = args.outdir+'/'
    allfiles = os.listdir(args.dir)
    libs = []
    for each in allfiles:
        if each.endswith('.fa'): libs.append(each)
        else: pass
    if args.outdir == None: args.outdir = args.dir
    else: pass
    try: os.mkdir(args.outdir)
    except OSError: pass

    #Make directory of all the info                                                                                                                
    loclib = {}
    for each in libs:
        library = each.split('.')[0]
        libfile = open(args.dir+each, 'r')
        for lines in libfile:
            info = lines.strip().split()
            if info[0].startswith('>'): locus = info[0][1:]
            else:
                seq = info[0]
                if locus not in loclib: loclib[locus] = {library:seq}
                else: loclib[locus][library] = seq
    for loc in loclib:
        locfile = open(args.outdir+loc+'.fa', 'w')
        for lib in loclib[loc]:
            locfile.write('>%s\n%s\n' % (lib, loclib[loc][lib]))
        locfile.close()


################
#  CALL FIXER  #
################

elif args.command == 'callFixer':
    if args.out.endswith('/'): pass
    else: args.out = args.out+'/'
    try: os.mkdir(args.out)
    except OSError: pass

    libs = []
    allFiles = os.listdir(args.bam)
    for each in allFiles:
        if each.endswith('.sorted.bam'): libs.append(each.split('.')[0])
        else: pass
    for lib in libs:
        print 'Making pileup file for library '+lib
        os.system('samtools mpileup -f %s%s.fa %s%s.sorted.bam > %s%s.pileup' % (args.ref, lib, args.bam, lib, args.out, lib))
        os.system('rm %s%s.fa.fai' % (args.ref, lib))
        print 'Calling assembly fixer for library '+lib
        os.system('perl %s FIX -p %s%s.pileup -a %s%s.fa -f %s%s -i 0.05'% (args.script, args.out, lib, args.ref, lib, args.out, lib))
        os.system('rm %s%s.pileup' % (args.out, lib))
        print 'Library %s done!' % lib
