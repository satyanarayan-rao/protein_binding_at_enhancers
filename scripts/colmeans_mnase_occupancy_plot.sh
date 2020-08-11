#!/bin/bash
# $1: input tsv file
# $2: column_id - 3, for example
# $3: output pdf 
# $4: output png 
# $5: input base gplt
# $6: output gplt 
# $7: plot title
# $8: bed slop value 

column_ids=`awk '{print $v}' v=$2 $1 | sort | uniq` 
plt_file=
plt_command="plot "
to_remove=""

for c in $column_ids
do
    grep -w $c"$" $1 > ${1}.${c}.tmp
	  plt_command=`echo "$plt_command '${1}.${c}.tmp' using 1:2 with lines title '$c' lw 2, "`
	  to_remove=`echo "$to_remove ${1}.$c.tmp " `
done

cp $5 $6
cat << EOM  >> $6 
set output '$3'
set title noenhanced '$7'
set xlab 'Distance from peak center[bp]'
set ylab 'E.O.M'
set xtics ('-${8}' -${8}, '0' 0, '${8}' ${8} )  
$plt_command 
EOM

gnuplot $6

convert -density 150 -background white -alpha remove $3 $4
# cleanup 
rm $to_remove 
