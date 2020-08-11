#!/bin/bash 
# $1: input file - peak center bed follwed by fragment bed
# $2: cluster id
# $3: left flank 
# $4: right flank
# $5: output file

awk '{if ( (($2 - $8) >= lflank) && (($9 - $2) >= rflank) ) {print $0}}' lflank=$3 rflank=$4 $1 | grep "\^${2}"  > ${5}
