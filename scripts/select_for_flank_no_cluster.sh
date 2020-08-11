#!/bin/bash
# $1: input file - peak center bed follwed by fragment bed
# $2: left flank
# $3: right flank
# $4: output file
awk '{if ( (($2 - $8) >= lflank) && (($9 - $2) >= rflank) ) {print $0}}' lflank=$2 rflank=$3 $1 > ${4}

