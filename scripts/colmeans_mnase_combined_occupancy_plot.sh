#!/bin/bash
# $1: input short tsv (comb_50)
# $2: input long tsv  (comb_147)
# $3: input cluster cnt tsv
# $4: input base gplt
# $5: output gnuplt file
# $6: output pdf
# $7: output png
# $8: clusted id  to plot  
# $9: plot title 

grep -w "Cluster${8}" $1 > ${1}.${8}.tsv 
grep -w "Cluster${8}" $2 > ${2}.${8}.tsv 
n=`awk '{if ($1==v) {print $2}}' v=$8 $3` 

cp $4 $5 

cat << EOM >> $5
set output  '$6'
set title noenhanced '$9 (n = $n)'  
set y2tics
set ytics nomirror
set y2tics textcolor rgb 'red'
set xlab 'Distance from peak center [bp]'
set ylab 'E.O.M' 
plot '${1}.${8}.tsv' using 1:2 w l title 'Short (50bp)' lw 2 lc rgb 'blue', '${2}.${8}.tsv' using 1:2 w l title 'Long (147bp)' lw 2 lc rgb 'red' axis x1y2
EOM

gnuplot  $5 
convert -density 150 -background white -alpha remove $6 $7
# cleanup
rm ${1}.${8}.tsv
rm ${2}.${8}.tsv
