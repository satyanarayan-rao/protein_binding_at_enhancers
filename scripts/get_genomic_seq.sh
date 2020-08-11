#!/bin/bash
# $1: overlapping bed gz 
# $2: output bed file gz 0-coordinate based
# $3: output bed file gz 1-coordinate based

zcat $1 | awk '{print $1"\t"$2-1"\t"$3"\t"$4}' | gzip - > ${2}.query.bed.gz # subtracting 1 to makw it zero coordinate ; please see `fasta_to_alignment_and_methylation_vector.md` 

bedtools getfasta -fi ~/data/ucsc/dm/dm3/dm3.fa  -bed  ${2}.query.bed.gz -tab -name+ | awk -F '[-:\t]' '{print $3"\t"$4"\t"$5"\t"$1"\t.\t"$NF}' | gzip - > ${2}
#bedtools getfasta -fi ~/data/ucsc/dm/dm3/dm3.fa  -bed  ${2}.query.bed.gz -tab -name+ | awk -F '[-:\t]' '{print $3"\t"$4+1"\t"$5"\t"$1"\t.\t"$NF}' | gzip - > ${3}

zcat $2 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}'  | gzip - > $3 

