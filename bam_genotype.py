#!/usr/bin/env python3
"""
Call genotypes using freebayes.

eg: bam_genotype.py -r /data/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa -t targets.gff /data/seq/AgamP4_coluzzii/04SELI0001.bam

"""
import sys
import os
import time
import argparse
import collections
import configparser
import shlex
import subprocess
import tempfile
import re
import codecs

import vcf


def num2str(val, none_val='0'):
    return none_val if val is None else str(val)

### Main ######################################################################

def Main(argv=None):
#    print("START", file=sys.stderr)
    start_tic = time.time()

    # parse cfg_file argument
    conf_parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter,
                                        add_help=False) # turn off help so later parse (with all opts) handles it
    conf_parser.add_argument('-c', '--cfg-file', type=argparse.FileType('r'), help='(optional) config file specifiying options/parameters')
    args, remaining_argv = conf_parser.parse_known_args(argv)

    # read the config as the argument defaults if given
    if args.cfg_file:
        cfg = configparser.ConfigParser()
        cfg.read_string("[DEFAULTS]\n"+args.cfg_file.read())
        defaults = dict(cfg.items("DEFAULTS"))
        defaults['cfg_file'] = args.cfg_file
        # special handling of paratmeters that need it like lists
        if( 'bam_files' in defaults ): # bam_files needs to be a list
            defaults['bam_files'] = [ x for x in defaults['bam_files'].split('\n') if x and x.strip() and not x.strip()[0] in ['#',';'] ]
    else:
        defaults = {}

    # Parse rest of arguments
    aparser = argparse.ArgumentParser(description=__doc__, parents=[conf_parser], formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    aparser.add_argument('bam_files', metavar='bam-file', nargs='*')
#    aparser.add_argument('-d', '--min-depth', type=int,
#                help="Min depth required per-sample to make a genotype call", default=2)
#    aparser.add_argument('-D', '--min-hom-depth', type=int,
#                help="Min per-sample depth required to call if all reads/observations are the same. if <= min-depth, min-depth is used.", default=5)
    aparser.add_argument('-P', '--single-population', action='store_true',
                help="All samples are assumed to come from a single population.  Default behaviour is to analyze each sample as its own population unless a --populations file is provided.")
    aparser.add_argument('-t','--targets',
                required= 'targets' not in defaults,
                help="File of target sites to call.  Either a GFF file, or else a VCF/BED-like (1-based) file with 'chrom position ...' on each line.")
    aparser.add_argument('-p','--populations_file',
                help="Populations file passed to freebayes.  Each line is a sample and a population which it is part of.")
    aparser.add_argument('-r','--reference',
                required= 'reference' not in defaults,
                help="Reference FASTA file (shoudl be same used for making/mapping bam_files)")
    aparser.add_argument('--ploidy', type=int, default='2',
                help="ploidy value to pass to freebayes")
    aparser.add_argument('-A','--cnv-map',
                help="copy-number map filename to pass to freebayes")
    aparser.add_argument('--freebayes-exec',
                help="Freebayes executable file",
                default='freebayes')
    aparser.add_argument('--samtools-exec',
                help="Samtools executable file",
                default='samtools')
    aparser.add_argument('--freebayes-options', help="Options used when running freebayes",
                default='--standard-filters --use-mapping-quality --max-complex-gap -1 --report-monomorphic')
    aparser.add_argument('-0','--zero-based-targets', help="Targets file is 0-based (like BED) as oppsoed to 1-based (like VCF)", action='store_true', default=False)
#    aparser.add_argument('--missing-char', help="Character to indicate missing/no-call", default='.')
#    aparser.add_argument('-V', '--output-vcf', action='store_true', help="Just output raw VCF results from freebayes")
    aparser.add_argument('-K', '--keep-tempfiles', action='store_true', help="keep temporary files")
#    aparser.add_argument('-x', '--no-partial-calls', action='store_true', help="Don't output partial calls (like A/-)")
    aparser.add_argument('-v', '--verbose', action='count', help="Increase verbosity level")
    aparser.set_defaults(**defaults)

    # process options/arguments
    args = aparser.parse_args(remaining_argv)

    args.output_vcf = True # @TCC TEMP OVERRIDE

    # @TCC TEMP -- print out all the options/arguments
    for k,v in vars(args).items() :
        print(k,":",v, file=sys.stderr)
#    sys.exit(1)

    # read the targets gff file and write out simplified targets file for freebayes
    targets_by_pos = ReadTargetsFile(args.targets, args.zero_based_targets)
    targets_tempfilename = WriteTargetsTempfile(targets_by_pos, args.keep_tempfiles)

    # get the sample names from the bam files
    samples = []
    for bam_file in args.bam_files:
        proc = subprocess.Popen([args.samtools_exec, 'view', '-H', bam_file], stdout=subprocess.PIPE)
        pat = re.compile(r'@RG\s.*\sSM:([^\s]+)')
        for line in iter(proc.stdout.readline,b''):
            m = pat.match(line.decode('utf-8'))
            if m :
                samples.append(m.group(1))

    # handle populations file or lack thereof
    if( args.populations_file ):
        pops_filename = args.populations_file
    else:
        if( args.single_population ):
            pops_filename = None
        else:
            # make a populations file with each sample in it's own pop
            pops_fh = tempfile.NamedTemporaryFile(mode='w', dir='./', delete=not args.keep_tempfiles)
            for i, samp in enumerate(samples):
                print(samp, "pop{0:04d}".format(i), file=pops_fh)
            pops_fh.flush()
            pops_filename = pops_fh.name

    # Run freebayes
    cmd = []
    cmd.append(args.freebayes_exec)
    cmd.extend(shlex.split(args.freebayes_options))
    cmd.extend(['-f', args.reference])
    cmd.extend(['-t', targets_tempfilename])
    cmd.extend(['--ploidy', str(args.ploidy)])
    if( args.cnv_map ):
        cmd.extend(['-A', args.cnv_map])
    if( pops_filename ):
        cmd.extend(['--populations', pops_filename])
    for bam_file in args.bam_files:
        assert os.path.isfile(bam_file)
        cmd.extend(['-b', bam_file])
    print("Running freebayes:", file=sys.stderr)
    print(" ".join(cmd), file=sys.stderr)

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)

    # Just output the vcf
    if( args.output_vcf ):
        for line in iter(proc.stdout.readline,''):
            sys.stdout.write(line)
        sys.exit(0)

    ## Parse the vcf to make calls (@TCC Not working and should probably be made separate)
    vcf_reader = vcf.Reader(iter(proc.stdout.readline,''))
    data = {}
    alleles = {}
    for record in vcf_reader:
        pos = "{0:s}:{1:010d}".format(record.CHROM,record.start)
        alleles[pos] = ','.join([str(x) for x in record.alleles])
        for samp in vcf_reader.samples:
            key = samp+' '+pos
            # make a list of the observations/read-counts.  [ref,alt1,alt2,...]
            obs = [ record.genotype(samp)['RO'] ]
            if( record.genotype(samp)['AO'] is None or isinstance(record.genotype(samp)['AO'], (int, float)) ):
                obs.append(record.genotype(samp)['AO'])
            else: # alternate observations is a list
                obs.extend(record.genotype(samp)['AO'])
            obs = [ 0 if x is None else x for x in obs ] # replace None with 0
            obs_str = ','.join([num2str(x) for x in obs])
            depth = sum(obs)
            #print(obs)

            gt_bases = record.genotype(samp).gt_bases
            if( gt_bases is not None ):
                gt_bases = gt_bases.split('/')
            #qual = ','.join([num2str(record.genotype(samp)['QR']), num2str(record.genotype(samp)['QA'])])
            # additional filters go here:
            if( depth < args.min_depth ):
                gt_bases = None
            if( gt_bases is not None and args.min_hom_depth > args.min_depth ):
                if( all(x == gt_bases[0] for x in gt_bases) and depth < args.min_hom_depth ):
                    if( args.no_partial_calls ):
                        gt_bases = None
                    else:
                        for i in range(1,len(gt_bases)):
                            gt_bases[i] = args.missing_char
            if( gt_bases is None ):
                gt_bases_str = args.missing_char
            else:
                gt_bases_str = '/'.join(gt_bases)
            #print(key, obs, depth, gt_bases, file=sys.stderr)
            data[key] = [ obs_str, gt_bases_str ]

    # output header lines
    # ...variant labels
    tmp = ['#']
    for pos,t in targets_by_pos.items():
        tmp.append(t['ID'])
        tmp.append(re.sub('\s+', '_', "{chrom}:{start:010d};{info}".format(**t)))
    # ...alleles
    tmp = ['#']
    for pos,t in targets_by_pos.items():
        tmp.append('-')
        tmp.append(alleles[pos])
    print('\t'.join(tmp))
    # ...column lables
#    tmp = ['#']
#    for pos,t in targets_by_pos.items():
#        tmp.append('RO,AO')
#        tmp.append('GT')
#    print('\t'.join(tmp))

    # output data
    for samp in vcf_reader.samples:
        tmp = [samp]
        for pos in targets_by_pos:
            key = samp+' '+pos
            if not key in data:
                tmp.extend([str(None),str(None)])
            else:
                tmp.extend(data[key])
        print('\t'.join(tmp))

    # Cleanup and end normally
    #os.unlink(targets_tempfilename) # remove the targets tempfile
    print("END", "%0.2f"%(time.time()-start_tic),"secs elapsed", file=sys.stderr)


def ReadTargetsFile(filename, zero_based_targets=False):
    filetype = 'vcf'
    _, ext = os.path.splitext(filename)
    if( ext.lower() == '.gff' ):
        filetype = 'gff'
    else:
        with open(filename,'r') as fh:
            line = fh.readline()
            if bool(re.match('##gff-version', line, re.I)) :
                filetype = 'gff'
    if( filetype == 'gff' ):
        return ReadGFFFile(filename, zero_based_targets)
    else:
        return ReadVCFLikeTargetFile(filename, zero_based_targets)


def ReadGFFFile(filename, zero_based_input=False):
    """Note: GFF start positions are 1-based and end is inclusive
        output will be 0-based
    """
    targets_by_pos = collections.OrderedDict()
    with open(filename,'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#') :
                continue
            f = line.split('\t')
            tmp = dict([ #collections.OrderedDict([
                    ('chrom', f[0]),
                    ('source', f[1]),
                    ('feature', f[2]),
                    ('start', int(f[3])),
                    ('end', int(f[4])),
                    ('length', int(f[4])-int(f[3])+1),
                    ('score', f[5]),
                    ('strand', f[6]),
                    ('frame', f[7]),
                    ('info', f[-1]),
                    ])
            if( not zero_based_input ): # convert from 1-based?
                tmp['start'] -= 1
                tmp['end'] -= 1
            # make sure start comes before end
            if tmp['start'] > tmp['end']:
                _ = tmp['start']
                tmp['start'] = tmp['end']
                tmp['end'] = _
            # ID?
            m = re.search(r'ID=([^;]*);',tmp['info'])
            if( m ):
                tmp['ID'] = m.group(0)[3:-1]
            else:
                tmp['ID'] = None
                print('WARNING: No ID for GFF entry:', file=sys.stderr)
                print('\t',line, file=sys.stderr)
            targets_by_pos["{0:s}:{1:010d}".format(tmp['chrom'],tmp['start'])] = tmp
    return targets_by_pos


def ReadVCFLikeTargetFile(filename, zero_based_input=False):
    """output will be 0-based"""
    targets_by_pos = collections.OrderedDict()
    with open(filename,'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#') :
                continue
            if line.startswith('browser') :
                continue
            if line.startswith('track') :
                continue
            f = line.split('\t')
            if( '=' in f[0] ):
                continue
            tmp = dict([ #collections.OrderedDict([
                    ('chrom', f[0]),
                    ('start', int(f[1])),
                    ('end', int(f[1])),
                    ('ID', "{0:s}:{1:010d}".format(f[0],int(f[1]))),
                    ('info', '-'),
                    ])
            if( not zero_based_input ): # convert from 1-based?
                tmp['start'] -= 1
            targets_by_pos["{0:s}:{1:010d}".format(tmp['chrom'],tmp['start'])] = tmp
    return targets_by_pos


def WriteTargetsTempfile(targets_by_pos, keep_tempfiles=False):
    # write the targets to a temp file for freebayes to use
    global targets_fh # so it is only deleted when program exits
    targets_fh = tempfile.NamedTemporaryFile(mode='w', dir='./', delete=not keep_tempfiles)
    #targets_fh = sys.stdout
    first_line = True
    for target_pos in sorted(targets_by_pos):
        t = targets_by_pos[target_pos]
        if( first_line ):
            print('\t'.join([ t['chrom'], '1', '2', 'dummy' ]), file=targets_fh) # dummy line needed due to freebayes bug
            first_line = False
        print('\t'.join([ t['chrom'], str(t['start']), str(t['end']+1), str(t['ID']) ]), file=targets_fh)
    targets_fh.flush()
    return targets_fh.name



#########################################################################
# Main loop hook... if run as script run main, else this is just a module
if __name__ == "__main__":
    sys.exit(Main(argv=None))
