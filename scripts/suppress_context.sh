#!/bin/bash
bam_file=$1
suppressed_bam=$2
samtools view -H ${bam_file} > ${bam_file}.tmp.header
samtools view ${bam_file} | python scripts/suppress_hch_and_dgd_contexts.py | cat ${bam_file}.tmp.header - | samtools view -Sb -o ${suppressed_bam} - 
