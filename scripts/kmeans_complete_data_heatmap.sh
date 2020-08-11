#!/bin/bash
# $1: input file
# $2: base gnuplt command
# $3: output gnuplt file
# $4: eps
# $5: pdf 


python scripts/kmeans_complete_data_heatmap.py $1 $2 $3 $4 "$6" $7

gnuplot $3

ps2pdf $4 $5

