#!/bin/bash 
bam_file=$1
flag=$2
tot_lines=$3

samtools view  $bam_file | awk '{if($2==flag){print $0}}' flag=$flag | awk '{print $10"\t"$(NF-2)}' | sed 's/XM:Z://g' | head -n $tot_lines
