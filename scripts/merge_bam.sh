#!/bin/bash
# $1: input list of bams "a.bam b.bam c.bam"
# $2: merged bam file
rm -f tmp/merge_list.${3}
for i in `echo $1`
do 
    echo $i >>  tmp/merge_list.${3}
done
bamtools merge -list tmp/merge_list.${3} -out $2
