#!/bin/bash
# $1: per_orange_hist_png 
# $2: hist_of_all_lengths_cap_200_combined_png 
convert +append $1 $2  ${1}.83.tmp.png
convert +append $3 $4 ${3}.83.tmp.png
convert +append $8 $9 ${8}.83.tmp.png
convert -density 150 -background white -alpha remove ${5} ${5%.pdf}.png 
convert -density 150 -background white -alpha remove ${6} ${6%.pdf}.png 

convert +append ${5%.pdf}.png ${6%.pdf}.png ${5%.pdf}.83.png 
convert -append ${5%.pdf}.83.png  ${3}.83.tmp.png ${1}.83.tmp.png ${8}.83.tmp.png  $7


convert +append ${1} ${3} ${8} ${7}.right.top.png 
convert +append ${5%.pdf}.png ${7}.right.top.png ${7}.top.png

convert +append ${2} ${4} ${9} ${7}.right.bottom.png 
convert +append ${6%.pdf}.png ${7}.right.bottom.png ${7}.bottom.png
# cleanup
rm ${5%.pdf}.83.png  ${1}.83.tmp.png ${3}.83.tmp.png ${8}.83.tmp.png
