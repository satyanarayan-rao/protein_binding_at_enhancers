import os 
import sys
from collections import defaultdict
import re
bed_fp = open(sys.argv[1]) 
loc_to_peak_map_fp = open(sys.argv[2])
out_fp = open(sys.argv[3], "w")
# bed_fp: footprint_min_size_selection/suppressed_merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_lf_15_rf_15_tf_footprint_min_length_10_wobble_gap_1.bed
# chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147	.	..........................................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF........................................................	..X..Z.....Z.............H..............................................H...Z.......H....H............................z....z...x....x....z....z.....z....h.....................................x........z...H.......z...................H...................Z......H.......z....................	GGCTGCGTTTGCGTGGAAAGAAGAGCAAAAATATTTTAATTTTAATTGGAAAGGGAATGGGAATTGGAATTGCTTTCGGGGGGGCATTGCAAATTTTATATTTTTAATTATTTTAATTTGAAATGAGTAGTGTTGTGTGGTGTGGTTTTGGTGTTAATAAATCAAATTTATTGTATTTATATTTAGAAGAGTTTTATGTGTGTGCATTTTTTTGATTAGTTTGTTTTTTTTGCTTTTTTTCTAGATTTTAGACGGGGGGCCTTTTTTTGTTTTTTTTTTTTTATTATT


# loc_to_peak_map_fp : input_bed/reordered_20clust_147_50PNE_open_v1_peak.bed
# chr3R	19975479	19975979	19975740	8	-	peak_4746	38.90	SingleZ	4
peak_read_cnt_dict = defaultdict(lambda : defaultdict(lambda : 0 ))
loc_to_peak_map_dict = {} 
loc_to_strand_dict = {} 
for line in loc_to_peak_map_fp: 
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    loc_id = line[0:d_loc[0]] + ":" + line[d_loc[2] + 1: d_loc[3]] + "-" + line[d_loc[2] + 1: d_loc[3]] 
    peak_id = line[d_loc[5] + 1: d_loc[6]]
    loc_to_peak_map_dict[loc_id] = peak_id
    loc_to_strand_dict[loc_id] = line[d_loc[4] + 1: d_loc[5]]

for line in bed_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    chr_loc_and_cl = line[d_loc[2] + 1: d_loc[3]] 
    sam_flag = line[d_loc[8] + 1: d_loc[9]].split("`")[-1]  
    peak_read_cnt_dict[chr_loc_and_cl][sam_flag]  += 1 

header = "%20s\t%40s\t%10s\t%10s\t%5s\t%3s"%("peak_id", "chr_loc", "99/147", "83/163", "mnase_cl", "strand")
out_fp.write (header + "\n")
for k in peak_read_cnt_dict:
    mnase_cl = k.split ("^")[-1]
    chr_loc = k.split ("^")[0]
    peak_id = loc_to_peak_map_dict[chr_loc]
    strand  = loc_to_strand_dict[chr_loc]
    to_write = "%20s\t%40s\t%10s\t%10s\t%5s\t%3s"%(peak_id,chr_loc, peak_read_cnt_dict[k]["99~147"],
                             peak_read_cnt_dict[k]["83~163"], mnase_cl, strand)
    out_fp.write(to_write + "\n")
