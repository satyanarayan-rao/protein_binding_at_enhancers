#!/bin/bash 
# $1: bam file
# $2: bedgz file
samtools view $1 | python scripts/prepapre_methylation_string.py | gzip - > $2 
