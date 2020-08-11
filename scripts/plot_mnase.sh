#!/bin/bash
# $1: input raw signal
# $2: input base raw gnuplt file
# $3: input eom signal
# $4: input base eom gnuplt file
# $5: input extact signal 
# $6: output raw eps
# $7: output raw gplt
# $8: output eom eps 
# $9: output eom gplt 
# ${10}: params flank bp annotation
# ${11}: raw pdf
# ${12}: eom pdf

# copy base gnuplot files to output gnuplot files 
cp $2 $7
cp $4 $9  

# Total number of lines in the matrix (apply to both raw and eom) 
tot_lines=`wc -l $5 | awk '{print $1}'`
#get total number of columns 
tot_cols=`zless $1 | head -1 | awk '{print NF}'` 
tot_cols=`expr $tot_cols - 1` # subtract 1 to get exact columns (1st one is bed details) 
half=`echo "$tot_cols / 2" | bc`
# Get coordinates to plot lines on heatmap 
lines_to_draw=`zcat $1 | awk '{print $1}' | awk -F'^' '{print $NF}' | uniq -c | awk '{sum+=$1; print sum}' | sed '$d'`
cat <<EOT >> $7
set output '$6'
set xra [1:$tot_cols] 
set yra [1:$tot_lines] 
set xtics ("-$half" 1, "0" $half, "+$half" $tot_cols)
EOT
cnt=1
for line_loc in $lines_to_draw
do
cat <<EOT >> $7
set arrow $cnt from 1, $line_loc to $tot_cols, $line_loc nohead lc rgb "#000000" front
EOT
cnt=`expr $cnt + 1`
done 
cat <<EOT >> $7
set xlab "Distance from enhancer center (bp)"
set cbra [0:5]
set cblabel "Raw MNase signal"
set title "DM enhancer sites (n = $tot_lines)"
plot '<zcat $1' matrix with image notitle
EOT
cat <<EOT >> $9
set output '$8'
set xra [1:$tot_cols] 
set yra [1:$tot_lines] 
set xlab "Distance from enhancer center (bp)"
set xtics ("-$half" 1, "0" $half, "+$half" $tot_cols)
EOT
cnt=1
for line_loc in $lines_to_draw
do
cat <<EOT >> $9
set arrow $cnt from 1, $line_loc to $tot_cols, $line_loc nohead lc rgb "#000000" front
EOT
cnt=`expr $cnt + 1`
done 
cat <<EOT >> $9
set cbra [0:5]
set cblabel "Enrichment over mean (MNase signal)"
set title "DM enhancer sites (n = $tot_lines)"
plot '<zcat $3' matrix with image notitle
EOT

gnuplot $7 
gnuplot $9 

ps2pdf $6 ${11}
ps2pdf $8 ${12}

