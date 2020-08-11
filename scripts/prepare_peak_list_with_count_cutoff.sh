#!/bin/bash 

sam_99_cutoff=$2
sam_83_cutoff=$3

awk 'NR>1{if ($3>=sam_99_cutoff && $4>=sam_83_cutoff){print $1"\t"$5}}' sam_99_cutoff=${sam_99_cutoff} sam_83_cutoff=${sam_83_cutoff} $1  
