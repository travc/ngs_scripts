#!/bin/bash

REF=$1

#echo "Trim blank lines (in place) and index '$REF'?"
#read -p "Continue (y/N)?"
#[ "$(echo $REPLY | tr [:upper:] [:lower:])" == "y" ] || exit

# trim blank lines
read -p "Trim blank lines of $REF [in place] (y/N)?"
if [ "$(echo $REPLY | tr [:upper:] [:lower:])" == "y" ] ; then
    echo
    echo "Trimming blank lines from $REF"
    #sed -i '/^$/d' $REF
fi

echo
echo "Picard tools CreateSequenceDictionary..."
java -jar /usr/local/picard-tools/CreateSequenceDictionary.jar R=$REF O=${REF%.*}.dict

echo
echo "Samtools faidx..."
samtools faidx $REF

echo
echo "BWA index..."
bwa index $REF

echo "Making BLAST db (formatdb)..."
formatdb -p F -o T -i $REF
