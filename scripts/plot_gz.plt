reset
set terminal postscript enhanced color size 6, 6 font 'Arial, 15'
set output 'check.eps' 
set palette defined ( 0 "#e0ecf4", 0.5 "#9ebcda", 1 "#8856a7")
#set palette defined ( 0 "#8856a7", 0.5 "#9ebcda", 1 "#e0ecf4")
set xra [1:501]
set yra [1:2731]
set cbra [30:100]
plot '<zcat ./tandem_matrix/S2_tandem_R1_2_CpG_and_S2_tandem_R1_2_GpC_to_50PNE_open_v1_chr_added_tandem_matrix.csv.gz'  matrix with image notitle
