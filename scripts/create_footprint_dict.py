import os
import sys
import pickle
import re
# head -1 footprint.bed
#chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147`99~147	.	........F...................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF......................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.....FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.....FFFFFFFFFFFFFFFFFFFFFFF.....FFF.....FFFFFFF......................................

# head -1 merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_7_lf_22_rf_19_methylation_matrix_for_site_peak_3484_nclust_3_sam_flag_83~163.tsv 
# chr3R:5416975-5416975^7`SRR3133326.727267_727267/1_overlapping`83~163`83~163#4	...................FFFFFFFFFFFFFF.....FFFF

input_fp = open(sys.argv[1]) 
out_fp = open(sys.argv[2], "wb")
footprint_dict = {}
for line in input_fp:
    # example line:
    # chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147	.	..........................................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF........................................................	..X..Z.....Z.............H..............................................H...Z.......H....H............................z....z...x....x....z....z.....z....h.....................................x........z...H.......z...................H...................Z......H.......z....................	GGCTGCGTTTGCGTGGAAAGAAGAGCAAAAATATTTTAATTTTAATTGGAAAGGGAATGGGAATTGGAATTGCTTTCGGGGGGGCATTGCAAATTTTATATTTTTAATTATTTTAATTTGAAATGAGTAGTGTTGTGTGGTGTGGTTTTGGTGTTAATAAATCAAATTTATTGTATTTATATTTAGAAGAGTTTTATGTGTGTGCATTTTTTTGATTAGTTTGTTTTTTTTGCTTTTTTTCTAGATTTTAGACGGGGGGCCTTTTTTTGTTTTTTTTTTTTTATTATT
    tmp = {}   
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    key = line[d_loc[2] + 1: d_loc[3]] + "`" + line[d_loc[8] + 1:d_loc[9]]
    tmp["start"] = line[d_loc[6]+1:d_loc[7]] 
    tmp["end"] = line[d_loc[7]+1:d_loc[8]] 
    tmp["footprint"] = line[d_loc[10] + 1: d_loc[11]]  
    tmp["mvec"] = line[d_loc[11] + 1: d_loc[12]] 
    tmp["bs_seq"] = line[d_loc[12] + 1: len(line) - 1] 
    footprint_dict[key] = tmp

pickle.dump (footprint_dict, out_fp) 
out_fp.close()
input_fp.close()
