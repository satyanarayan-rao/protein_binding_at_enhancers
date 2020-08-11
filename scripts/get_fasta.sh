#!/bin/bash
# $1: input bed file
# $2: genome fasta file
# $3: fasta output file 
bedtools getfasta -fi $2 -bed $1 -fo $3 -name 
