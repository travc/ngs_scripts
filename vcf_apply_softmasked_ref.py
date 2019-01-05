#!/usr/bin/env python3
"""Filter out lines from a vcf file based where the fasta reference is lower-case or 'N'"""

import sys
import os
import time
import argparse

import pyfaidx


def exit(retcode=0):
    # cleanup and exit
    global g_start_tic
    print("END {:0.2f} secs elapsed".format(time.time()-g_start_tic), file=sys.stderr)
    sys.exit(retcode)


def Main(argv=None):
    global g_start_tic
    g_start_tic = time.time()

    aparser = argparse.ArgumentParser(description=__doc__)
    aparser.add_argument('-r', '--ref', type=argparse.FileType('r'), help="Reference fasta (presumably softmasked); default is to infer from VCF header", default=None)
    aparser.add_argument('-v', '--vcf', type=argparse.FileType('r'), help="Input VCF file to filter", default='-')
    args = aparser.parse_args()

    infh = args.vcf

    if args.ref is None:
        # read through the header until we get the reference
        ref_filename = None
        for line in infh:
            if line.startswith('##reference='):
                ref_filename = line.strip().split('=')[1]
                sys.stdout.write(line)
                break
            if line.startswith('##'):
                sys.stdout.write(line)
            else:
                print("ERROR: Could not find '##reference=' line in header", file=sys.stderr)
                exit(1)
        args.ref = open(ref_filename,'r')

    print("REFERENCE = '{}'".format(args.ref.name), file=sys.stderr)

    # load the reference
    ref_fasta = pyfaidx.Fasta(args.ref.name)

    # read through rest of vcf
    missing_chroms = []
    for line in infh:
        if line.startswith('#'):
            sys.stdout.write(line)
        else:
            ch, pos = line.split()[0:2]
            try:
                base = str(ref_fasta[ch][int(pos)-1])
            except KeyError as e:
                if not ch in missing_chroms:
                    print("WARNING: {}:{} not found in ref.  Suppressing further warnings for this chrom".format(ch, pos), file=sys.stderr)
                    missing_chroms.append(ch)
                    continue
            if base.isupper() and base != 'N':
                sys.stdout.write(line)
            else:
                pass # skip
                #print(str(ref_fasta[ch][int(pos)-1]))
                #sys.stdout.write(line)
                #break

    exit(0)




#########################################################################
# Main loop hook... if run as script run main, else this is just a module
if __name__ == "__main__":
    sys.exit(Main(argv=sys.argv))
