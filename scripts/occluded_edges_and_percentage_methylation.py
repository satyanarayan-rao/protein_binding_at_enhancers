import sys
import pickle
import re 
from collections import defaultdict
from capital_peracentage_and_unbound_stretches import capital_percentage_and_stretch_of_unbound

inp_fp = open(sys.argv[1])

out_fp_complete_mvec_stats = open(sys.argv[2], "w") 
out_fp_complete_mvec_stats_pkl = open(sys.argv[3], "wb") 
occluded_dna_dict = defaultdict(dict) 

for line in inp_fp:
# example line:
# chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147	.	..........................................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF........................................................	..X..Z.....Z.............H..............................................H...Z.......H....H............................z....z...x....x....z....z.....z....h.....................................x........z...H.......z...................H...................Z......H.......z....................	GGCTGCGTTTGCGTGGAAAGAAGAGCAAAAATATTTTAATTTTAATTGGAAAGGGAATGGGAATTGGAATTGCTTTCGGGGGGGCATTGCAAATTTTATATTTTTAATTATTTTAATTTGAAATGAGTAGTGTTGTGTGGTGTGGTTTTGGTGTTAATAAATCAAATTTATTGTATTTATATTTAGAAGAGTTTTATGTGTGTGCATTTTTTTGATTAGTTTGTTTTTTTTGCTTTTTTTCTAGATTTTAGACGGGGGGCCTTTTTTTGTTTTTTTTTTTTTATTATT

    d_loc = [m.start() for m in re.finditer("\t",line)]
    m_vec = line[d_loc[11] + 1: d_loc[12]] 
    # make the key the as: "chr2L:11806729-11806729^20`SRR3133326.952549_952549/1_overlapping`99~147"
    k = line[d_loc[2] + 1: d_loc[3]] + "`" + line[d_loc[8] + 1:d_loc[9]] 
    # Implementing the occlulded logic 
    # getting the left and right occluded edges
    mvec_len, total_on_complete_read, pcap, occluded_dna_stretches =\
       capital_percentage_and_stretch_of_unbound (m_vec)
    max_of_edges = max(occluded_dna_stretches[0], occluded_dna_stretches[2])
    to_write = "\t".join([k, str(mvec_len), str(total_on_complete_read),
                str(pcap), "\t".join(map(str,occluded_dna_stretches)),
                str(max_of_edges)]) 
    out_fp_complete_mvec_stats.write(to_write + "\n") 
    occluded_dna_dict[k]["total_vec_len"] = mvec_len
    occluded_dna_dict[k]["total_letters_on_read"] = total_on_complete_read
    occluded_dna_dict[k]["percentage_capital"] = pcap
    occluded_dna_dict[k]["max_of_edge"] = max_of_edges 
    if occluded_dna_stretches[0] >= occluded_dna_stretches[2]:
        occluded_dna_dict[k]["start"] = 0
        occluded_dna_dict[k]["end"] = occluded_dna_stretches[0] - 1 
    else:
        occluded_dna_dict[k]["start"] = mvec_len - occluded_dna_stretches[2]
        occluded_dna_dict[k]["end"] = mvec_len - 1 
    
     
pickle.dump(occluded_dna_dict, out_fp_complete_mvec_stats_pkl) 
out_fp_complete_mvec_stats.close() 
out_fp_complete_mvec_stats_pkl.close()

 
