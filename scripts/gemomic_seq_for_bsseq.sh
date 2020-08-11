#!/bin/bash
# $1: bsseq  file
# $2: output genomic file
awk '{print $7"\t"$8-1"\t"$9"\t"$10}' $1 |  bedtools getfasta -fi ~/data/ucsc/dm/dm3/dm3.fa -bed - -tab -name+ | awk -F '[-:\t]' '{print $3"\t"$4+1"\t"$5"\t"$1"\t.\t"$NF}' > ${2}.tmp.bed

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $1 > ${1}.enh
paste -d"\t" ${1}.enh ${2}.tmp.bed  > ${2}

