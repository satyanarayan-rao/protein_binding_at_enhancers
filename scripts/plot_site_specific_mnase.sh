#!/bin/bash
# $1: short fragment e.csv.gz mnase input
# $2: long fragment e.csv.gz mnase input
# $3: base gnuplt file
# $4: output plot tsv file
# $5: output gplt file
# $6: output pdf file
# $7: string for site to search
# $8: span from the center
# $9: left vertical line
# ${10}: right vertical line
# ${11}: plot title

cp ${3} $5 

#echo $1
#echo $2
#echo $3
#echo $4
#echo $5
#echo $6
#echo $7
#echo $8
#echo $9
#echo ${10}
#echo ${11}

zless $1 | grep "$7" | tr '\t' '\n' | awk 'NR>1{print $0}' > ${4}.s.tmp
tot_lines=`wc -l ${4}.s.tmp | awk '{print $1}'`
half=`echo "${tot_lines}/2" | bc`
left_from_center=`echo 0 - ${10} | bc`
awk '{print NR - half"\t"$0}' half=$half ${4}.s.tmp  | awk '{if ($1 >=left_from_center && $1 <=right_from_center){print $0}}' left_from_center=-$9 right_from_center=${10} > ${4}.short.tmp 
max_of_short=`sort -k2,2n ${4}.short.tmp | tail -1 | awk '{print $2}'`
zless $2 | grep "$7" | tr '\t' '\n' | awk 'NR>1{print $0}' > ${4}.l.tmp

tot_lines=`wc -l ${4}.l.tmp | awk '{print $1}'`
half=`echo "${tot_lines}/2" | bc`
awk '{print NR - half"\t"$0}' half=$half ${4}.l.tmp  | awk '{if ($1 >=left_from_center && $1 <=right_from_center) {print $0}}' left_from_center=-$9 right_from_center=${10} > ${4}.long.tmp 
max_of_long=`sort -k2,2n ${4}.long.tmp | tail -1 | awk '{print $2}'` 
max_of_both=`echo -e "${max_of_short}\t${max_of_long}" | awk '{if($1>=$2){print $1}else{print $2}}'`

paste ${4}.short.tmp ${4}.long.tmp | awk '{print $1"\t"$2"\t"$4}' > $4 

cat <<EOT >> $5
set output "$6" 
set title "${11}"
set ylab "E.O.M" 
set y2range [0:${max_of_long}]
set yrange[0:${max_of_short}]
set y2tics textcolor rgb "red"
set ytics nomirror
set style arrow 1 nohead ls 4 lw 4  lc rgb "black"
set arrow 1 arrowstyle 1 from -${9}, 0 to -${9}, $max_of_both
set arrow 2 arrowstyle 1 from ${10}, 0 to ${10}, $max_of_both 

plot "$4"  using 1:2 w l title "Short fragments (<50)" lc rgb "blue" lw 4 axis x1y1, "$4" using 1:3 w l title "Nuc Fragments (147)" lc rgb "red" lw 4 axis x1y2
EOT

gnuplot $5 
Rscript scripts/mnase_plot_for_top_panel.R $4 ${12}

# cleanup 
#rm ${4}.s.tmp ${4}.l.tmp ${4}.short.tmp ${4}.long.tmp
