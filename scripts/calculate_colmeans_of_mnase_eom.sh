#!/bin/bash
# $1: eom csv gz 
# $2: nclust
# $3: slop
# $3: output file
# zcat mnase_based_on_occupancy/sm.sub.comb_50_c_mapped_on_suppressed_merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_ocm_selected_vplot_cluster_3_peaks_with_lf_15_rf_15_condition_cond4.kmeans_nclust_3.merged.slop_500.csv.gz | grep "\`1" | perl -lane 'for $c (1..$#F){$t[$c] += $F[$c]}; END{for $c (1..$#t){print $t[$c]/$.}}' | cat -n  |
to_sub=`echo $3 + 1 | bc`
rm ${4}.tmp
for cl_id in `seq $2`
do
    zcat $1 | grep "\`${cl_id}" | perl -lane 'for $c (1..$#F){$t[$c] += $F[$c]}; END{for $c (1..$#t){print $t[$c]/$.}}' | cat -n  | awk '{print $1 - to_sub"\t"$2"\tCluster"cl_id }' to_sub=${to_sub} cl_id=${cl_id}  >> ${4}.tmp
done 
mv ${4}.tmp $4 

