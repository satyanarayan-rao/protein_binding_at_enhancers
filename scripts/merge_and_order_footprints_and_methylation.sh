#!/bin/bash

# $1 : footprint 99 
# $2 : footprint 83 
# $3 : methylation 99
# $4 : methylation 83
# $5 : footprint merged and ordered
# $6 : methylation merged and ordered


python scripts/sam_flag_specific_merge.py $1 $2 ${5}.tmp.tsv
python scripts/sam_flag_specific_merge.py $3 $4 ${6}.tmp.tsv

#cat $1 $2 > ${5}.tmp.tsv
#cat $3 $4 > ${6}.tmp.tsv
#
#python scripts/order_footprints_by_length_and_direction_left_ordered.py ${5}.tmp.tsv ${5}.tmp2.tsv
#python scripts/order_footprints_by_length_and_direction_left_ordered.py ${6}.tmp.tsv ${6}.tmp2.tsv
python scripts/footprint_dot_to_digit_vec.py ${5}.tmp.tsv ${5} 
python scripts/methylation_dot_to_digit_vec.py ${6}.tmp.tsv ${6} 

# cleanup
rm ${5}.tmp.tsv ${6}.tmp.tsv 
