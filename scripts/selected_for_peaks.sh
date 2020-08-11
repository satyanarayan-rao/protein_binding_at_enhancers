#!/bin/bash 
# $1: kmeans tsv 
# $2: input bed 
# $3: output bed

# example lines for $1
#file_name_without_sam_flag	total_reads	total_naked	total_tf	total_nuc	occup_naked	occup_tf	occup_nuc	cl_id
#binding_labeled_reads_for_peak/suppressed_merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_6_lf_15_rf_15_site_peak_2057_binding_label.cond4.tsv	152	104	9	39	68.42	5.92	25.66	6


# example lines for $2 
#chr3R	19975479	19975979	19975740	8	-	peak_4746	38.90	SingleZ	4
#chrX	3636439	3636939	3636843	7	-	peak_3219	15.50	SingleZ	4
rm ${3}.${4}.tmp.bed
awk 'NR>1 {print $1"\t"$NF}' $1 | awk -F"_site_"  '{print $NF}' | awk -F '[_\t]' '{print $1"_"$2"\t"$NF}' > ${1}.${4}.tmp.tsv 

while read peak_id cl_id 
do
   grep -w "${peak_id}" $2 | awk '{print $1"\t"$4"\t"$4+1"\t"peak"`"cl_id"\t.\t"$6}' peak=${peak_id} cl_id=${cl_id} >>  ${3}.${4}.tmp.bed
done < ${1}.${4}.tmp.tsv 

bedtools slop -b $4 -g $5 -i ${3}.${4}.tmp.bed > $3
# cleanup 
rm ${1}.${4}.tmp.tsv
