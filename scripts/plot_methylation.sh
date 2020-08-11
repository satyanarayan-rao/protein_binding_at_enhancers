#!/bin/bash
# $1: input tandem matrix
# $2: base gnuplt file
# $3: output eps
# $4: output pdf 
# $5: output gnuplt file

cp $2 $5
tot_lines=`zcat $1 | wc -l`
tot_cols=`zless $1 | head -1 | awk '{print NF}'`
tot_cols=`expr $tot_cols - 1`
half=`echo "$tot_cols / 2" | bc`

cat <<EOT >> $5
set output '$3'
set xra [1:$tot_cols] 
set yra [1:$tot_lines]
set cbra [0:20] 
set xlab "Distance from enhancer center (bp)"
set cblabel "CpG/GpC methylation percentage" 
set title "{/Arial-Italic dm} enhancer sites (n = $tot_lines)"
set xtics ("-$half" 1, "0" $half, "+$half" $tot_cols)
plot '<zcat $1' matrix with image notitle
EOT

gnuplot $5 

ps2pdf $3 $4 
