#!/bin/bash

# $1 : annotated footprint 
# $2 : non_annotated verbose file (have more )
python scripts/order_footprints_by_length_and_direction_left_ordered.py $1 ${3}
awk 'NR%3 == 2' $2 > ${4}.tmp
python scripts/annotate_cobinding_using_footprint.py $1  $4
