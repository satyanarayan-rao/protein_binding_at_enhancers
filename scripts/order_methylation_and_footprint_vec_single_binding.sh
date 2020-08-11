#!/bin/bash

# $1 : verbos file with edges as E
# $2 : output footprint ordered_footprint vec
# $3 : output ordered methylation vector


awk 'NR%3==1' $1 > ${2}.fp.tsv 
awk 'NR%3==2' $1 > ${3}.mvec.tsv 


python scripts/order_footprints_by_length_and_direction_left_ordered.py ${2}.fp.tsv ${2}
python scripts/order_footprints_by_length_and_direction_left_ordered.py ${3}.mvec.tsv ${3}

# cleanup

rm ${2}.fp.tsv ${3}.mvec.tsv
