#!/bin/bash
# $1: intersected bed file
# $2: input bam file (superset)
# $3: output bam file (subset)
# $4: sample name
# $5: bed name
# $6: picard jar path

awk '{print $4}' $1 | awk -F'/' '{print $1"/"$2}' > tmp/tmp.${4}.${5}.read_id_list
java -jar $6 FilterSamReads I=$2 O=$3 READ_LIST_FILE=tmp/tmp.${4}.${5}.read_id_list FILTER=includeReadList
# cleanup
#rm tmp/tmp.${4}.${5}.read_id_list  

