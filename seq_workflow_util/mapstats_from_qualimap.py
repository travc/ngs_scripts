#!/usr/bin/env python3
"""
"""
import sys
import os
import time
import argparse
import glob
import collections
import configparser
import itertools
import shlex
import pipes
import subprocess
import tempfile
import re
import errno
import fileinput
from collections import OrderedDict

import numpy

### Constants #################################################################
DEFAULT_CONFIG = None #'map_stats.cfg' # Default config filename; None for none; Useful for testing

EXCLUDE_DIRS = ['qualimap'] # dirs to exclude from sample dir glob matching

QUALIMAP_WGS_OUTDIR = 'qualimap'

# a bit of python 2 to 3 compatibility
try:
  basestring
except NameError:
  basestring = str

################################################################################

#### Utility Functions (could be moved to a module) ####

def path_split_all(path):
    tmp = path
    split_path = []
    while( tmp ):
        tmp, tmp2 = os.path.split(tmp)
        if( tmp2 != '' ):
            split_path.insert(0, tmp2)
        if( tmp == os.sep ):
            split_path.insert(0, tmp)
            break
    return split_path

## math-ish ##

def hist_to_cumulative(inhist, normalize=True, reverse=False):
        # covert a histogram to a cumdist
        tsum = numpy.zeros_like(inhist[0,1:])
        out = numpy.zeros_like(inhist)
        for i,t in enumerate(inhist):
            tsum += t[1:]
            out[i,0] = t[0]
            out[i,1:] = tsum
        if( normalize ):
            out[:,1:] = numpy.divide(out[:,1:],tsum)
            tsum = numpy.ones_like(tsum) # normalized sum is all ones
        if( reverse ):
            out[:,1:] = numpy.subtract(tsum,out[:,1:])
        return out


## ConfigParser helper stuff (doesn't normally need any changing)

def print_args(args, out_file_handle=sys.stderr, prefix=''):
    """print out the args (Namespace) created by argparse
       'prefix' is just added to beginning of each line"""
    for k,v in vars(args).items() :
        print(prefix+str(k),":",str(v), file=out_file_handle)

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
## end: ConfigParser helper stuff


def transpose_data_dict(t):
    """converts a dict with 'fileds' and 'vals' (where vals are lists of lists)
    and coverts it into a dict where the keys are the fields and the values are the associate vals"""
    foo = {}
    for i,row in enumerate(t['vals']):
        for j,val in enumerate(row):
            field = t['fields'][j]
            if( not field in foo ):
                foo[field] = []
            foo[field].append(val)
    return foo

################################################################################

## YL: fix this part to parse qualimap/genome_results.txt ###
def s2i(s):
    return int(s.replace(',','').strip().split()[0])
def s2f(s):
    return float(s.replace(',','').strip().split()[0])

def parse_qualimap(filename, contig_list):
    out = OrderedDict()
    with open(filename) as fh:
        for line in fh:
            #l = line.strip().split(maxsplit=3)
            l = line.strip().split(sep='=')

            if l[0].startswith('number of reads'):
                out['total_reads'] = s2i(l[1]) # NOTE: counts reads with secondary mapping twice!
            elif l[0].startswith('number of mapped reads'):
                out['mapped_reads'] = s2i(l[1])
            elif l[0].startswith('median insert size'):
                out['insert_size'] = s2i(l[1])
            elif l[0].startswith('mean mapping quality'):
                out['mean_mapping_quality'] = s2f(l[1])
            elif l[0].startswith('GC percentage'):
                out['GC_content'] = float(l[1].replace('%','').strip())
            elif l[0].startswith('mean coverageData'):
                out['WG_cov'] = float(l[1].replace('X','').strip())
            elif l[0].startswith('std coverageData'):
                out['WG_cov_std'] = float(l[1].replace('X','').strip())
            ## regions
            else:
                for k in contig_list:
                    if l[0].startswith(k+'\t'):
                        l2 = l[0].split('\t')
                        out[k+'_cov'] = float(l2[3])
                        out[k+'_cov_std'] = float(l2[4])

        #out['reads'] = out['total_reads']-out['secondary_reads']
    return out

################################################################################


### Main ######################################################################
def exit(retcode=0):
    # cleanup and exit
    global g_start_tic
    print("END", "%0.2f"%(time.time()-g_start_tic),"secs elapsed", file=sys.stderr)
    sys.exit(retcode)


def Main(argv=None):
#    print >>sys.stderr, "START"
    global g_start_tic
    g_start_tic = time.time()

    # parse cfg_file argument
    conf_parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=CustomArgparseHelpFormatter,
                                        add_help=False) # turn off help so later parse (with all opts) handles it
    conf_parser.add_argument('-c', '--cfg-file', type=argparse.FileType('r'),
            help="Config file specifiying options/parameters.\nAny long option can be set by remove the leading '--' and replace '-' with '_'",
            default=DEFAULT_CONFIG)
    args, remaining_argv = conf_parser.parse_known_args(argv)

    # build the config (read config files)
    if args.cfg_file:
        cfg = configparser.SafeConfigParser()
        cfg.optionxform=str # make the ConfigParser options case sensitive
        cfg_lines = itertools.chain(("[DEFAULTS]",), args.cfg_file)
        cfg.read_file(cfg_lines)
        defaults = dict(cfg.items("DEFAULTS"))
        # special handling of paratmeters that need it like lists
        for k in ['sample_dirs', 'contigs']:
            if( k in defaults ):
                defaults[k] = [ x for x in defaults[k].split('\n') if x and x.strip() and not x.strip()[0] in ['#',';'] ]
    else:
        defaults = {}

    # Parse rest of arguments with a new ArgumentParser
    aparser = argparse.ArgumentParser(description=__doc__, parents=[conf_parser], formatter_class=CustomArgparseHelpFormatter)

    aparser.add_argument('sample_dirs', metavar='sample-dir', nargs='*', help='sample directory')
    aparser.add_argument('--input-bam', default='merged_markdup.bam',
                help='name of the bam file in each sample_dir to analyze (needed for idxstats)')
    aparser.add_argument('--gfc-thresholds', default='75,50,25',
                help='comma separated list of % values to report coverage where genome fraction is >=')
    aparser.add_argument('--contigs', default='all',
                help='list (comma separated) of contigs to report coverage data for. "all" or "none" for all contigs or none respectively')
    #aparser.add_argument('-G', '--graphic-outfile', default=None,
    #            help='Filename of graphical report to generate. "none" or "" for no graphical report.')
    #aparser.add_argument('-C', '--csv-outfile', default=None,
    #            help='Filename of csv report to generate. "none" or "" for no graphical report.')

    aparser.add_argument('-v', '--verbose', action='count', help="Increase verbosity level")

    aparser.set_defaults(**defaults)

    # process options/arguments
    args = aparser.parse_args(remaining_argv)

    # custom/more complex argument parsing errors
    if not hasattr(args, 'verbose') or args.verbose is None:
        args.verbose = '0'
    # numbers which are given as strings so stuff like '1e5' works
    for arg_name in ['verbose', ]:
        if( getattr(args, arg_name) is not None ):
            if( getattr(args, arg_name).lower() == 'none' ):
                setattr(args, arg_name, None)
            else:
                setattr(args, arg_name, int(round(float(getattr(args, arg_name)))))
    ## expand potential globs (wildcards) in a file list
    # also checks that file exists
    sample_dirs = []
    for f in args.sample_dirs:
        tmp = glob.glob(f)
        tmp = [x for x in tmp if os.path.basename(os.path.normpath(x)) not in EXCLUDE_DIRS]
        if not tmp:
            continue
        if( not tmp ):
            raise IOError('Cannot find dir(s) matching "'+str(f)+'"')
        sample_dirs.extend(tmp)

    # comma separated list args
    args.gfc_thresholds = [float(x) for x in args.gfc_thresholds.split(',')]

    if( args.verbose > 0 ):
        print_args(args)

    # contigs may need to be determined from data if 'all' is given
    if isinstance(args.contigs, basestring):
        args.contigs = args.contigs.split(',')
    if( args.contigs[0] == 'none' ):
        contig_list = None
    elif( args.contigs[0] == 'all' ):
        contig_list = []
    else:
        contig_list = args.contigs
    ##
    data = [] # results by sample

    #### ####
    header_done = False
    for sample_dir in sorted(sample_dirs):

        # do we need to fill out the contig_list?
        if contig_list == []:
            ## from qualimap "Coverage per contig" output
            #tmp = glob.glob(sample_dir+'/qualimap/genome_results.txt')
            #contig_list = []


            #for tmp2 in tmp:
                #contig_list.append(re.search('.*flagstat_(.*)\.out', tmp2).group(1))
            #contig_list.sort()
            contig_list = ['1','2','3','MT']  # for Aedes aegypti
            print('contig_list =',contig_list, file=sys.stderr)

        # get sample name from the full path
        tmp = sample_dir
        while True:
            tmp, sample_name = os.path.split(tmp)
            if( sample_name != '' or not tmp ):
                break

        print('*'*78,"\nProcessing:", sample_name, sample_dir, file=sys.stderr)


        # start adding info for this sample to data
        data.append(OrderedDict())
        d = data[-1]
        d['name'] = sample_name
        d['path'] = sample_dir
        d.update(parse_qualimap(os.path.join(sample_dir, 'qualimap', 'genome_results.txt'), contig_list))
        #d['trim'] = read_trim_stats(sample_dir)

        for region in contig_list:
            d[region] = OrderedDict()

        #
        #for k,v in d.items():
        #    print(k,v, file=sys.stderr)

        ## header if first sample ##
        if not header_done:
            header_done = True
            print(  "#name",
                    "num raw reads", # actually reads survivied trimming
                    "num mapped reads",
                    "% mapped",
                    "med insert size",
                    "WG med cov",
                    "WG med cov std",
                    "raw reads/med_cov",
                    sep='\t',
                    end='\t',
                    file=sys.stdout)
            for region in contig_list:
                print(
                    region+" cov",
                    region+" cov std",
                    sep='\t',
                    end='\t',
                    file=sys.stdout)
            sys.stdout.write('\n')

        ## Sample summary output ##
        print(  d['name'],
                d['total_reads'], 
                d['mapped_reads'], 
                "{:4.2f}%".format(100*(d['mapped_reads'])/float(d['total_reads'])), # mapped % (excluding secondary)
                int(d['insert_size']), # median insert size
                #"{:.1f}".format(d['insert_size']['is_mean']), # mean insert size
                int(d['WG_cov']), # median insert size
                int(d['WG_cov_std']), # median insert size
                sep='\t',
                end='\t',
                file=sys.stdout)
        for region in contig_list:
            print(
                #"{:4.2f}%".format(100*d[region]['frac_covered']),
                #"#{:.0f}".format(d['trim']['num_raw_reads']/d[region]['cov_50%']),
                "{:.2f}".format(d[region+'_cov']),
                "{:.2f}".format(d[region+'_cov_std']),
                sep='\t',
                end='\t',
                file=sys.stdout)
        sys.stdout.write('\n')

        ####################################################################################################


    # Cleanup and end normally
    exit()


#########################################################################
# Main loop hook... if run as script run main, else this is just a module
if __name__ == "__main__":
    sys.exit(Main(argv=None))
