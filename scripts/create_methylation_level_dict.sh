#!/bin/bash
# $1: bam file
# $2: output file
# $3: wildcard 
samtools view $1 | awk '{v=(NR-1)%2+1;print $1"/"v"\t"$3"\t"$4"\t"$14}'  | awk -F'XM:Z:' '{print $1""$2}' > tmp/tmp.${3}.tsv
python scripts/methylation_vector_to_dict.py tmp/tmp.${3}.tsv $2
# clean up
rm tmp/tmp.${3}.tsv
