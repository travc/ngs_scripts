#!/usr/bin/python3
"""convert a vcf to a sample-per-row tsv table.
simple implementation only really intended to work on SNPs.
"""

import sys
import os
import errno
import time
import argparse
import logging
from collections import OrderedDict
import vcf

MAX_LOGGING_LEVEL = logging.CRITICAL
DEFAULT_LOGGING_LEVEL = logging.INFO

def setup_logger(verbose_level):
    fmt=('%(levelname)s %(asctime)s [%(module)s:%(lineno)s %(funcName)s] :: '
            '%(message)s')
    logging.basicConfig(format=fmt, level=max((0, min((MAX_LOGGING_LEVEL,
                        DEFAULT_LOGGING_LEVEL-(verbose_level*10))))))

def Main(argv=None):
    tic_total = time.time()

    #parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', nargs='?', metavar='VCF', type=argparse.FileType('r'),
                        default='-',
                        help="VCF file to convert.  '-' for stdin (default).")
    parser.add_argument('-d', '--min-depth', type=int, default=0,
                        help="filter calls made with less than this FMT/DP (high-quality reads)")
    parser.add_argument('-c', '--report-counts', action='store_true',
                        help="Output supporting read counts")
    parser.add_argument('-F', '--soft-filter', action='store_true',
                        help="instead of dropping calls failing filter, just prefix '?' to them")
    parser.add_argument('-R', '--no-report-ref', action='store_true',
                        help="Do not report the reference as first data row (and ALT as second row).")

    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="increase logging verbosity")
    parser.add_argument('-q', '--quiet', action='count', default=0,
                        help="decrease logging verbosity")

    args = parser.parse_args()
    setup_logger(verbose_level=args.verbose-args.quiet)

    logging.info("vcf = "+str(args.vcf))
    logging.info("min_depth = "+str(args.min_depth))

    vcf_reader = vcf.Reader(args.vcf)
    samples = vcf_reader.samples

    dat = []
    pos = []
    ref_call = []
    alt_call = []
    for record in vcf_reader:
        logging.info(str(record))
        pos.append(str(record.CHROM)+':'+str(record.POS))
        d = []
        ref_call.append(record.REF)
        alt_call.append(','.join((str(x) for x in record.ALT)))
        for sample in samples:
            c = record.genotype(sample).gt_bases
            if ( c is None or
                 record.genotype(sample)['DP'] < args.min_depth ):
                if args.soft_filter:
                    #c = '?'+''.join(sorted(str(c).split('/')))
                    c = '?'+str(c)
                else:
                    c = '?'
                #c = '?DP='+str(record.genotype(sample)['DP'])
            else:
                c = str(c)
                #c = ''.join(sorted(str(c).split('/')))
            if args.report_counts:
                c += '\t'
                c += str(record.genotype(sample)['RO'])
                c += ','
                #print(len(record.ALT), str(record.genotype(sample)['AO']))
                try:
                    c += ','.join((str(_) for _ in record.genotype(sample)['AO']))
                except TypeError:
                    c += str(record.genotype(sample)['AO'])
            d.append(c)
        dat.append(d)

    print('# file:', args.vcf.name)
    print('# min_depth:', args.min_depth)

    if args.report_counts:
        print('', *[_2 for _1 in zip(pos,['RO,AO']*len(pos)) for _2 in _1], sep='\t')
        if not args.no_report_ref:
            print('REF', end='\t')
            print(*ref_call, sep='\t\t')
            print('ALT', end='\t')
            print(*alt_call, sep='\t\t')
        for i,s in enumerate(samples):
            print(s, *[dat[j][i] for j in range(len(pos))], sep='\t')

    else:
        print('', *pos, sep='\t')
        if not args.no_report_ref:
            print('REF', *ref_call, sep='\t')
            print('ALT', *alt_call, sep='\t')
        for i,s in enumerate(samples):
            print(s, *[dat[j][i] for j in range(len(pos))], sep='\t')

    logging.info("Done: {:.2f} sec elapsed".format(time.time()-tic_total))


#########################################################################
# Main loop hook... if run as script run main, else this is just a module
if __name__ == "__main__":
    sys.exit(Main(argv=None))

