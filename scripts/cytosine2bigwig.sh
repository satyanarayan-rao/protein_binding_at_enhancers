#!/bin/bash
# $1: input cov file
# $2: output bigwig file
gunzip -c $1 | awk '{print $1"\t"$2"\t"$3+1"\t"$4}' | sort -k1,1 -k2,2n - > ${1}.tmp
bedGraphToBigWig ${1}.tmp $3 $2

# cleanup

rm ${1}.tmp
