#!/usr/bin/env python
"""
counts indicies from illumina fastq file
"""
import sys
import os
import time
import argparse
import glob
import collections
import ConfigParser
import shlex
import pipes
import subprocess
import tempfile
import re
import errno
import fileinput
import gzip

### Constants #################################################################

################################################################################

def num2str(val, none_val='0'):
    return none_val if val is None else str(val)

class ConfigFakeSecHead(object):
    def __init__(self, fp, section='DEFAULTS'):
        self.fp = fp
        self.sechead = '['+str(section)+']\n'
    def readline(self):
        if self.sechead:
            try: return self.sechead
            finally: self.sechead = None
        else: return self.fp.readline()

class CustomArgparseHelpFormatter(argparse.HelpFormatter):
    """Help message formatter for argparse
        combining RawTextHelpFormatter and ArgumentDefaultsHelpFormatter
    """

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])

    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help


################################################################################


### Main ######################################################################
def exit(retcode=0):
    # cleanup and exit
    global g_start_tic
    print >>sys.stderr, "END", "%0.2f"%(time.time()-g_start_tic),"secs elapsed"
    sys.exit(retcode)


def Main(argv=None):
    print >>sys.stderr, "START"
    global g_start_tic
    g_start_tic = time.time()

    # parse cfg_file argument
    conf_parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=CustomArgparseHelpFormatter,
                                        add_help=False) # turn off help so later parse (with all opts) handles it
    conf_parser.add_argument('-c', '--cfg-file', type=argparse.FileType('r'), help="Config file specifiying options/parameters.\nAny long option can be set by remove the leading '--' and replace '-' with '_'")
    args, remaining_argv = conf_parser.parse_known_args(argv)

    # build the config (read config files)
    if args.cfg_file:
        cfg = ConfigParser.SafeConfigParser()
        cfg.readfp(ConfigFakeSecHead(args.cfg_file))
        defaults = dict(cfg.items("DEFAULTS"))
        # special handling of paratmeters that need it like lists
        if( 'fastq_files' in defaults ): # fastq_files needs to be a list
            defaults['fastq_files'] = [ x for x in defaults['fastq_files'].split('\n') if x and x.strip() and not x.strip()[0] in ['#',';'] ]
    else:
        defaults = {}

    # @TCC TEMP WHILE DEVELOPING
    #defaults['fastq_files'] = ['/data/archive/2014-07-21/140715_HS3B/Undetermined_indices/Sample_lane1/lane1_Undetermined_L001_R1_001.fastq.gz']
    #defaults['fastq_files'] = ['/data/archive/2014-07-21/140715_HS3B/Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz']

    # Parse rest of arguments with a new ArgumentParser
    aparser = argparse.ArgumentParser(description=__doc__, parents=[conf_parser], formatter_class=CustomArgparseHelpFormatter)

    aparser.add_argument('fastq_files', metavar='fastq-file', nargs='*')
    aparser.add_argument('-n', '--num-reads', default=None, help="Number of reads to look at, default=all")
    aparser.add_argument('-N', '--num-barcodes', default=10, type=int, help="Number of (most frequent) barcode sequences to report on; 0 for all ")
    aparser.add_argument('-v', '--verbose', action='count', help="Increase verbosity level")

    aparser.set_defaults(**defaults)

    # process options/arguments
    args = aparser.parse_args(remaining_argv)

    # custom/more complex argument parsing errors
    # need fastq_files list
    if( not args.fastq_files ):
        aparser.error("Must provide fastq files as arguments and/or in the CFG_FILE (fastq_files parameter)")

    ## expand potential globs (wildcards) in the bam_files list
    # also checks that bam files exist
    fastq_files = []
    for f in args.fastq_files:
        tmp = glob.glob(f)
        if( not tmp ):
            raise IOError('Cannot find fastq file(s) matching "'+str(f)+'"')
        fastq_files.extend(tmp)

    counts = {}
    for filename in fastq_files:
        print >>sys.stderr, "#"*180
        print >>sys.stderr, "Reading:", f

        if args.num_reads:
            tmp = int(float(args.num_reads))

        with gzip.open(filename) as f:
            in_block_count = 0 # counter used to divide up blocks
            for l in f:

                if( in_block_count == 0 ): # header line is first line in block
                    #sys.stdout.write(l)
                    idx = l.split(':')[-1][:-1]
                    if( idx in counts ):
                        counts[idx] += 1
                    else:
                        counts[idx] = 1

                in_block_count += 1
                if( in_block_count > 3 ): # blocks are 4 lines long
                    in_block_count = 0
                    # only look at part of file?
                    if args.num_reads:
                        tmp -= 1
                        if( tmp < 0 ):
                            break

    sorted_list = [(k,v) for v,k in sorted( [(v,k) for k,v in counts.items()],reverse=True) ]

    if args.num_barcodes > 0:
        to_report = sorted_list[0:args.num_barcodes]
    else:
        to_report = sorted_list
    for v in to_report:
        print v

    # Cleanup and end normally
    exit()


## functions #############################################################################




#########################################################################
# Main loop hook... if run as script run main, else this is just a module
if __name__ == "__main__":
    sys.exit(Main(argv=None))
