#!/usr/bin/env python
"""
counts indicies from illumina fastq file
WARNING... INCOMPLETE
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
import itertools
import Levenshtein

### Constants #################################################################

N7_adapters = [
        ['N701', 'TAAGGCGA'],
        ['N702', 'CGTACTAG'],
        ['N703', 'AGGCAGAA'],
        ['N704', 'TCCTGAGC'],
        ['N705', 'GGACTCCT'],
        ['N706', 'TAGGCATG'],
        ['N707', 'CTCTCTAC'],
        ['N708', 'CAGAGAGG'],
        ['N709', 'GCTACGCT'],
        ['N710', 'CGAGGCTG'],
        ['N711', 'AAGAGGCA'],
        ['N712', 'GTAGAGGA'],
        ]
N5_adapters = [
        ['N501', 'TAGATCGC'],
        ['N502', 'CTCTCTAT'],
        ['N503', 'TATCCTCT'],
        ['N504', 'AGAGTAGA'],
        ['N505', 'GTAAGGAG'],
        ['N506', 'ACTGCATA'],
        ['N507', 'AAGGAGTA'],
        ['N508', 'CTAAGCCT'],
        ]

################################################################################

#def hamming(str1, str2):
#    assert len(str1) == len(str2)
#    return sum(c1 != c2 for c1, c2 in itertools.izip(str1, str2))

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

def checkAdapterDistances(adapters, calcDistance, max_distance, out=sys.stderr):
    out_width = 6
    foo = [ [0]*len(adapters) for _ in range(len(adapters)) ]
    for i,a in enumerate(adapters):
        for j,b in enumerate(adapters):
            foo[i][j] = calcDistance(a[1],b[1])
    if( out is not None ): # some fancy-ish output
        print >>out, ' '.rjust(out_width),
        for a in adapters:
                print >>out, str(a[0]).rjust(out_width),
        print >>out
        for i,a in enumerate(foo):
            print >>out, str(adapters[i][0]).rjust(out_width),
            for b in a:
                print >>out, str(b).rjust(out_width),
            print >>out
    for i,a in enumerate(foo):
        for j,b in enumerate(a):
            if( i != j and b <= max_distance ):
                print >>sys.stderr, "ERROR: Adapters", adapters[i], "and", adapters[j], "only have distance", b, ", but max-distance is set to", max_distance
                sys.exit(0)


def exhuastiveMatch(idx_seq, adapters, calcDistance, max_distance, warning_fails=True, verbose=0):
    """exhaustively looks for best match
        return None if there is a tie"""
    match = None
    tie = False
    for i,a in enumerate(adapters):
        d = calcDistance(idx_seq, a[1])
        if( d <= max_distance ):
            if( match is None ):
                match = [i,d]
            else:
                print >>sys.stderr, "WARNING: Multiple matches for",idx_seq
                if( warning_fails ):
                    sys.exit(2)
                if( d == match[1] ):
                    tie = True
                if( d < match[1] ):
                    tie = False
                    match = [i,d]
        if( tie ):
            match = None
        if( verbose >= 3 ):
            print >>sys.stderr, i, idx_seq, a[1], d
    return match


### Main ######################################################################
def exit(retcode=0):
    # cleanup and exit
    global g_start_tic
    print >>sys.stderr, "END", "%0.2f"%(time.time()-g_start_tic),"secs elapsed"
    sys.exit(retcode)


def Main(argv=None):
    print >>sys.stderr, "START"
    print >>sys.stderr, "WARNING--- INCOMPLETE SCRIPT"
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
    defaults['fastq_files'] = ['/data/archive/2014-07-21/140715_HS3B/Undetermined_indices/Sample_lane2/lane2_Undetermined_L002_R1_001.fastq.gz']

    # Parse rest of arguments with a new ArgumentParser
    aparser = argparse.ArgumentParser(description=__doc__, parents=[conf_parser], formatter_class=CustomArgparseHelpFormatter)

    aparser.add_argument('fastq_files', metavar='fastq-file', nargs='*')
    aparser.add_argument('-H', '--hamming', action='store_true', help="Use Hamming distance instead of Levenshtein distance")
    aparser.add_argument('-d', '--max-distance', type=int, default=1, help="Maximum distance (errors) allowed to when assigning a read")
    aparser.add_argument('-W', '--warning_fails', action='store_true', help="Exit if there is a warning")
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

    # distance method
    if( args.hamming ):
        calcDistance = Levenshtein.hamming
    else:
        calcDistance = Levenshtein.distance

    # sanity checking the index adapters distances from each other
    print >>sys.stderr, "N7 index adapter distances:"
    checkAdapterDistances(N7_adapters, calcDistance, args.max_distance)
    print >>sys.stderr, "N5 index adapter distances:"
    checkAdapterDistances(N5_adapters, calcDistance, args.max_distance)


    counts = {}
    foo_counter = 0
    match_count = 0
    for filename in fastq_files:
        print >>sys.stderr, "#"*180
        print >>sys.stderr, "Reading:", f

        with gzip.open(filename) as f:
            tmp = 4e6
            in_block_count = 0 # counter used to divide up blocks
            for l in f:

                if( in_block_count == 0 ): # header line is first line in block
                    #sys.stdout.write(l)
                    idx = l.split(':')[-1][:-1]
                    # count the indices
                    if( idx in counts ):
                        counts[idx] += 1
                    else:
                        counts[idx] = 1

                    assert len(idx) == 16, 'Only setup for dual-indexing (N7+N5 length should be 16)'
                    N7 = idx[0:8]
                    N5 = idx[8:]

                    N7_match = exhuastiveMatch(N7, N7_adapters, calcDistance, args.max_distance, args.warning_fails, verbose=args.verbose)
                    N5_match = exhuastiveMatch(N5, N5_adapters, calcDistance, args.max_distance, args.warning_fails, verbose=args.verbose)

                    if( args.verbose > 1 ):
                        print >>sys.stderr, N7, N5, N7_match, N5_match

                    if( N7_match is not None and N5_match is not None ):
                        print >>sys.stderr, "*** MATCH ***", N7, N5, N7_match, N5_match
                        print >>sys.stderr, "             ", N7_adapters[N7_match[0]][1], N5_adapters[N5_match[0]][1], N7_adapters[N7_match[0]][0], N5_adapters[N5_match[0]][0]
                        match_count += 1


                    foo_counter += 1
                    if( foo_counter > 1e6 ):
                        print >>sys.stderr, match_count, "matches"
                        sys.exit(0)


                in_block_count += 1
                if( in_block_count > 3 ): # blocks are 4 lines long
                    in_block_count = 0
                # @TCC temp only look at part of file
                tmp -= 1
                if( tmp < 0 ):
                    break

    sorted_list = [(k,v) for v,k in sorted( [(v,k) for k,v in counts.items()],reverse=True) ]

    for v in sorted_list[0:10]:
        print v

    # Cleanup and end normally
    exit()


## functions #############################################################################




#########################################################################
# Main loop hook... if run as script run main, else this is just a module
if __name__ == "__main__":
    sys.exit(Main(argv=None))
