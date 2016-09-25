import sys
import os
import numpy as np
import h5py
import allel
import pandas
from collections import OrderedDict

__version__ = '0.1.0'

## Utility Functions

def str2int(s):
    """return the int value of string, handles strings like 1e6 too"""
    rv = None
    try:
        rv = int(s)
    except ValueError:
        rv = int(float(s))
    return rv

def str2range(s):
    chrom = None
    start = 1
    stop = None
    tmp = s.split(':')
    chrom = tmp[0]
    if len(tmp)>1:
        if '-' in tmp[1]:
            tmp = tmp[1].split('-')
        else:
            tmp = tmp[1].split('..')
        start = str2int(tmp[0])
        if len(tmp)>1:
            stop = str2int(tmp[1])
    return (chrom, start, stop)

def test_str2range1():
    assert str2range('2L:1-1e6') == ('2L', 1, 1000000)
def test_str2range2():
    assert str2range('2L:1e6') == ('2L', 1000000, None)
def test_str2range3():
    assert str2range('2L') == ('2L', 1, None)
def test_str2range4():
    # very long ints should be preserved perfectly
    assert str2range('2L:1234567891011121314151617') == ('2L', 1234567891011121314151617, None)

# Run tests (@TCC remove for script)
test_str2range1()
test_str2range2()
test_str2range3()
test_str2range4()


## Functions for loading callset data

def samplenames2idxs(callset, samplenames):
    """
    Lookup samplenames/id in callset 
    returns list of ids (orderd same as in callset) 
    and the list indexes of those ids in the callset
    If samplenames == None, return all"""
    tmp = allCallsetSampleIds(callset)
    if samplenames == None:
        ids = tmp
        idxs = list(range(len(tmp)))
    else:
        idxs = [ tmp.index(x) for x in samplenames ]
        ids = [ tmp[x] for x in idxs ]
    return ids, idxs


def allCallsetSampleIds(callset):
    """return all sample ids from callset
    Warning: only looks at first chrom entry (assumes same samples for all chroms)
    """
    return [ _.decode('utf-8') for _ in callset[list(callset.keys())[0]]['samples'] ]


def loadGenotypeArray(callset,
                       range_string,
                       sample_idxs,
                       MIN_FMTDP=0,
                       MAX_MISSING=None,
                       MAX_AF=0,
                       FILTER_SNP=False,
                       FILTER_BIALLELIC=False,
                       FILTER_NON_SYNONYMOUS_CODING=False,
                       verbose=0,
                      ):
    """
    load the callset for given region into GenotypeArray
    do basic / first pass filtering as we go
    return ( position array, genotype array, filter boolean array )
    """
#    print('##',range_string,'#'*60)
    ch, start, stop = str2range(range_string)
    ac = None

    # make SortedIndex of positions so we can quickly locate_range
#    print(ch, end='\t')
    pos = allel.SortedIndex(callset[ch]['variants']['POS'])
#    print(pos.shape[0])
    if pos.shape[0] == 0: # return empty if nothing for chrom
        return [],[],[]

    # create the slice
    try:
        sl = pos.locate_range(start,stop)
        #print(sl)
        pos = pos[sl]
    except KeyError:
        pos = []

    if len(pos) == 0: # no loci in slice
        g = []
        flt = []
    else:

        g = allel.GenotypeChunkedArray(callset[ch]['calldata']['genotype'])[sl].take(sample_idxs, axis=1)
        if verbose >= 1:
            print(range_string, g.shape, sep='\t')

        num_loci_in = g.shape[0]
        flt = np.ones(num_loci_in, dtype=bool)

        # filter SNP
        if FILTER_SNP:
            flt_snp = callset[ch]['variants']['TYPE'][sl] == b'snp'
            flt = flt & flt_snp
#            print('=',np.count_nonzero(flt), 'passing previous filters & SNP')

        # filter NON_SYNONYMOUS_CODING
        if FILTER_NON_SYNONYMOUS_CODING:
            flt_nonsyn = callset[ch]['variants']['EFF']['Effect'] == b'NON_SYNONYMOUS_CODING'
            flt = flt & flt_nonsyn
#            print('=',np.count_nonzero(flt), 'passing previous filters & NON_SYNONYMOUS_CODING')

        # filter genotypes on FMT:DP
        if MIN_FMTDP > 0:
            genoflt_FMTDP = callset[ch]['calldata']['DP'][sl].take(sample_idxs, axis=1) < MIN_FMTDP
            g[genoflt_FMTDP] = [-1,-1]
            tmp_num_calls = g.shape[0]*g.shape[1]
            tmp = np.count_nonzero(genoflt_FMTDP)
#            print('{} genotype call of {} ({:02.2f}%) fail FMT:DP filter'.format(
#                    tmp, tmp_num_calls, 100*tmp/float(tmp_num_calls)))

        # filter max_missing (genotype calls)
        if MAX_MISSING is not None and MAX_MISSING < len(sample_idxs):
            flt_max_missing = np.sum(g.is_missing(), axis=1) <= MAX_MISSING
            tmp = num_loci_in - np.count_nonzero(flt_max_missing)
#            print('{} loci of {} ({:02.2f}%) have > {} missing genotypes'.format(
#                    tmp, num_loci_in, 100*tmp/float(num_loci_in), MAX_MISSING))
            flt = flt & flt_max_missing
#            print('=',np.count_nonzero(flt), 'passing previous filters & max_missing')

        # fliter biallelic
        if FILTER_BIALLELIC:
            if ac is None:
                ac = g.count_alleles()
            flt_biallelic = ac.allelism() == 2
            flt = flt & flt_biallelic
#            print('=',np.count_nonzero(flt), 'passing previous filters & biallelic')

        if MAX_AF is not None and MAX_AF > 0:
            if ac is None:
                ac = g.count_alleles()
            flt_max_AF = ac.to_frequencies().max(axis=1) <= MAX_AF
            flt = flt & flt_max_AF
#            print('=',np.count_nonzero(flt), 'passing previous filters & max_AF')

        # apply filters
        g = g.compress(flt, axis=0)
        pos = pos.compress(flt, axis=0)
#        print('\t',pos.shape[0])

    return pos, g, flt


def loadMultipleGenotypeArrays(callset,
                                sample_idxs=None,
                                ranges=None,
                                MIN_FMTDP=0,
                                MAX_MISSING=None,
                                MAX_AF=0,
                                FILTER_SNP=False,
                                FILTER_BIALLELIC=False,
                                FILTER_NON_SYNONYMOUS_CODING=False,
                                verbose=0,
                                ):
    """load the callset for each region into a dict (by region) of GenotypeArrays
    """

    pos_dict = OrderedDict()
    g_dict = OrderedDict()
    flt_dict = OrderedDict()

    if sample_idxs is None:
        sample_idxs = list(range(len(allCallsetSampleIds(callset))))
    if ranges is None:
        ranges = list(callset.keys())

    for rngstr in ranges:
        pos_dict[rngstr], g_dict[rngstr], flt_dict[rngstr] = loadGenotypeArray(callset,
                                                                rngstr,
                                                                sample_idxs,
                                                                MIN_FMTDP,
                                                                MAX_MISSING,
                                                                MAX_AF,
                                                                FILTER_SNP,
                                                                FILTER_BIALLELIC,
                                                                FILTER_NON_SYNONYMOUS_CODING,
                                                                verbose,
                                                                )
    return pos_dict, g_dict, flt_dict


def loadMultipleIntoSingleGenotypeArray(callset,
                                sample_idxs=None,
                                ranges=None,
                                MIN_FMTDP=0,
                                MAX_MISSING=None,
                                MAX_AF=0,
                                FILTER_SNP=False,
                                FILTER_BIALLELIC=False,
                                FILTER_NON_SYNONYMOUS_CODING=False,
                                verbose=0,
                                ):
    """loads ranges from callset
    concating results into single GenotypeArray
    sample_idxs is None will load all in callset
    ranges is None will load all in callset
    """
    if sample_idxs is None:
        sample_idxs = list(range(len(allCallsetSampleIds(callset))))
    if ranges is None:
        ranges = list(callset.keys())
    _, g_dict, _ = loadMultipleGenotypeArrays(callset,
                                sample_idxs,
                                ranges,
                                MIN_FMTDP,
                                MAX_MISSING,
                                MAX_AF,
                                FILTER_SNP,
                                FILTER_BIALLELIC,
                                FILTER_NON_SYNONYMOUS_CODING,
                                verbose,
                                )
    chroms = [ch for ch in g_dict.keys() if len(g_dict[ch])>0]
    if len(chroms) == 0: # No data returned
        all_g = []
    else:
        all_g = g_dict[chroms[0]]
        all_g = all_g.vstack(*[g_dict[ch] for ch in chroms[1:]])
    return all_g


## Functions for metadata handling

def loadMetaFile(meta_fn,
                callset_all_sample_ids=None,
                sample_ids_to_drop=None,
                sample_column_idx=0,
                header_lines=0,
                sep='\t'):
    """reads a tsv/csv file into a pandas dataframe
    reorders to match order in callset_all_sample_ids
    adds 'idx' column
    """
    meta = pandas.read_csv(meta_fn, sep=sep, header=header_lines, index_col=sample_column_idx, comment='#')
    # Removing samples if requested
    if sample_ids_to_drop:
        for s in sample_ids_to_drop: # Drop one at a time so we can tolerate already missing samples
            try:
                meta.drop((s), inplace=True)
            except ValueError as e:
                print('WARNING:', e, file=sys.stderr)
    # reorder to match order in callset (callset_all_sample_ids)
    meta = meta.reindex([x for x in callset_all_sample_ids if x in meta.index.values])
    sample_ids = meta.index.values
    sample_idxs = [callset_all_sample_ids.index(x) for x in sample_ids]
    meta['idx'] = pandas.Series(sample_idxs, index=meta.index)
    return sample_ids, sample_idxs, meta

