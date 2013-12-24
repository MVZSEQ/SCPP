#!/usr/bin/home python

#######################################################################
#                                                                     #
#                    Pool ABySS Assemblies                            #
#                                                                     #
#    This script concatenates the different kmers and kmer            #
#    coverage assemblies into one file to be run with the             #
#    final assebmly script.                                           # 
#                                                                     #
#    Script written by: Joshua Penalba (joshua.penalba@anu.edu.au)    #
#    Written on: 20 Oct 2013       Last Modified: 20 Oct 2013         #
#                                                                     #
#######################################################################

import argparse
import sys
import os

pool = argparse.ArgumentParser(description='Pools all the various kmers and kmer coverages to one file per library for the final assembly step')

pool.add_argument('-a', dest='assemblies', help='directory containing all assemblies (formatted for ABYSS assembly done in the SCPP pipeline)')
pool.add_argument('-o', dest='out', help='output directory')
pool.add_argument('-i', dest='infofile', help='library info file if batch processing of all ')

if len(sys.argv)==1:
    pool.print_help()
    sys.exit(1)
args = pool.parse_args()

if args.infofile:
    args.homedir = '/'.join(args.infofile.split('/')[:-1])
    args.homedir = args.homedir+'/'

if args.assemblies: 
    if not args.assemblies.endswith('/'): args.assemblies = args.assemblies+'/'

if args.out.endswith('/'): pass
else: args.out = args.out+'/'

### START OF THE SCRIPT

#Creates a list of paths that will be ran
pathlist = []
projs = set()
if args.infofile:
    infofile = open(args.infofile, 'r')
    infofile.readline()
    for lines in infofile:
        info = lines.strip().split()
        if info[3] not in projs:
            projs.add(info[3])
            pathlist.append(args.homedir+info[3]+'/AbyssAssembly')
        else: pass
else: pathlist.append(args.assemblies)

#Makes the output directories
if args.out:
    try: os.mkdir(args.out)
    except OSError: pass
else:
    for path in pathlist:
        projpath = '/'.join(path.split('/')[:-1])
        try: os.mkdir(projpath+'/FinalAssembly')
        except OSError: pass

#Concatenate the files
for path in pathlist:
    projpath = '/'.join(path.split('/')[:-1])
    allfiles = os.listdir(path)
    libs = set()
    for each in allfiles:
        lib = each.split('_')[0]
        libs.add(lib)
    libs = list(libs)
    for library in libs:
        if args.out: os.system('cat %s%s_*-contigs.fa > %s%s.fa' % (path, library, args.out, library))
        else: os.system('cat %s/%s_*-contigs.fa > %s/FinalAssembly/%s.fa' % (path, library, projpath, library))
