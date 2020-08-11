#!/bin/bash

# $1 :  annotated footprint vector file 
# $2 : NOT annotated verbose file 
# $3 : ordered footprint vec
# $4 : ordered methylation vector

awk 'NR%3==2' $2 > ${2}.mvec.tsv 

python scripts/annotate_mvec_using_footprints.py $1 ${2}.mvec.tsv ${4}.tmp.tsv 

python scripts/order_footprints_by_length_and_direction_left_ordered.py ${1} ${3}
python scripts/order_footprints_by_length_and_direction_left_ordered.py ${4}.tmp.tsv ${4}

# cleanup

rm ${2}.mvec.tsv ${4}.tmp.tsv
