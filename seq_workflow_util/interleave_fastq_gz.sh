#!/bin/bash
#
# interleave a pair of fastq.gz files and append any additional fastq.gz files to the end
# output stream is uncompressed
# typical usage: interleave_fastq_gz.py trim_pe1.fastq.gz trim_pe2.fastq.gz trim_sing.fastq.gz ... | bwa mem -p ...

SEP=`echo -e '\x1E'`

RP1=$1
shift
RP2=$1
shift

# interleave the first pair
paste -d '$SEP' <(zcat ${RP1} | paste -d '$SEP' - - - -) \
      <(zcat ${RP2} | paste -d '$SEP' - - - -) \
  | tr '$SEP' '\n'

# just add all the rest to the end
while (($#)); do
    RS=$1
    zcat ${RS}
    shift
done

