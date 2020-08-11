#!/bin/bash
# $1 - u plot 99 
# $2 - u plot 83
# $3 - dot plot 99
# $4 - dot plot 83
# $5 - out png

convert -density 150 -background white -alpha remove ${3} ${3%.pdf}.u.png
convert -density 150 -background white -alpha remove ${4} ${4%.pdf}.u.png
convert +append $1 ${3%.pdf}.u.png ${4%.pdf}.u.png $2 $5

# cleanup 
rm ${3%.pdf}.u.png ${4%.pdf}.u.png
