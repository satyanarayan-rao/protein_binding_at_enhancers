#!/bin/bash 
# $1: 99 flag pdf
# $2: 83 flag pdf 
# $3: scatter plot png 99 
# $4: scatter plot png 83 
# $5: hist footprint png 99
# $6: hist footprint png 83
# $7: hist orange png 99
# $8: hist orange png 83
# $9: output png
convert -density 150 -background white -alpha remove $1 ${1}.tmp.png 
convert -density 150 -background white -alpha remove $2 ${2}.tmp.png 
# fist append vertical 99 
convert -append $3 $5 $7 ${3}.append.png 
convert -append $4 $6 $8 ${4}.append.png 

convert +append ${3}.append.png ${1}.tmp.png ${2}.tmp.png  ${4}.append.png $9

# cleanup
rm ${1}.tmp.png ${2}.tmp.png ${3}.append.png ${4}.append.png 
