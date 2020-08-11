import os
import sys
import re
import pickle
from collections import defaultdict
from collections import OrderedDict
def get_real_footprint_length (m_vec, m_vec_start, m_vec_stop, complete_vec): 
    """
    extract complete footprint length from the selected methylation vector
    exmaple: 
    complete vec: . . . . F F F F . . . F F F . . . . . F F F F F 
    index       : 0 1 2 3 4 5 6 7 8 9 1011121314151617181920212223
    flank vec:  :             F F . . . F F F . . . . . F
                              6 7 8 9 10111213141516171819 
    output: [4, 3, 5]
    """
    # first get all footprint lengths in m_vec
    flen_lengths = [len(a) for a in m_vec.split(".")] 
    if m_vec.startswith("F"):
        flen_lengths[0] = flen_lengths[0] +  len(complete_vec[0:m_vec_start].split(".")[-1])
    if m_vec.endswith("F"):
        flen_lengths[-1] = flen_lengths[-1] +  len(complete_vec[m_vec_stop:].split(".")[0])
    # exclude all zeros
    return_list = [a for a in flen_lengths if a!=0] 
    # prepare a footprint length vector - at each index it will tell what is the footprint size that index it associated with
    first_f = False
    cnt = 0
    loc_first = []
    gap = True
    for c in m_vec:
        if c == 'F':
            first_f = True
            if gap == True:
                loc_first.append (cnt)
                gap = False
        else:
            first_f = False
            gap = True
        cnt +=1
    out_vec = [0]*len(m_vec)
    for i in range (len(loc_first)):
        for j in range(loc_first[i], len(m_vec)):
            if m_vec[j] == "F":
                out_vec[j] = return_list[i]
            else:
                break
    return return_list, out_vec
        

input_fp = open(sys.argv[1])
#sam_flag_dict = pickle.load(open(sys.argv[2], "rb"))
lflank = int(sys.argv[2])
rflank = int(sys.argv[3])
output_fp = open(sys.argv[4], "w")
real_length_dict = open (sys.argv[5], "wb")
footprint_dict_fp = open(sys.argv[6], "wb")
strand_agnostic_footprint_dict_fp = open(sys.argv[7], "wb")
strand_agnostic_fp = open(sys.argv[8], "w") 

footprint_length_at_bp_res_fp = open(sys.argv[9], "wb")
footprint_length_at_bp_res_dict = {}
footprint_length_at_bp_res_tsv_fp = open(sys.argv[10], "w")
real_footprint_length_dict =  OrderedDict()
footprint_dict = {}
strand_agnostic_footprint_dict = {}
#dbg_file = open("tmp/check.tsv", "w")

for line in input_fp:
    # oooooooooo
    # 0123456789
    #     45
    # 
    # 15994238 - 50 - 15994009  
    # line example:  
    # chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147	.	..........................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF..............FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF....................................	..X..Z..x..Z.............H.........h...h.h..............................Hh..Z.......H....H........h....hh..h....hh....z....z...x....x....z....z.....z....h...h....H.............h.....h........x........z...H..h.hh.z..hx...x...h..h.x..H....h..H......h....Z......HH.h.hh.z...h.h.hhh.h.h.hh...	GGCTGCGTTTGCGTGGAAAGAAGAGCAAAAATATTTTAATTTTAATTGGAAAGGGAATGGGAATTGGAATTGCTTTCGGGGGGGCATTGCAAATTTTATATTTTTAATTATTTTAATTTGAAATGAGTAGTGTTGTGTGGTGTGGTTTTGGTGTTAATAAATCAAATTTATTGTATTTATATTTAGAAGAGTTTTATGTGTGTGCATTTTTTTGATTAGTTTGTTTTTTTTGCTTTTTTTCTAGATTTTAGACGGGGGGCCTTTTTTTGTTTTTTTTTTTTTATTATT

    delim_loc  = [m.start() for m in re.finditer("\t", line)]
    strand = line [delim_loc[4] + 1: delim_loc[5]]
    if strand == "-":
        lflank = int (sys.argv[3])
        rflank = int (sys.argv[2])
    elif strand == "+":
        lflank = int (sys.argv[2])
        rflank = int (sys.argv[3])
    m_vec_start = int(line[delim_loc[0]+1:delim_loc[1]]) - int(line[delim_loc[6]+1:delim_loc[7]]) -  lflank
    m_vec_stop = int(line[delim_loc[0]+1:delim_loc[1]]) - int(line[delim_loc[6]+1:delim_loc[7]]) + rflank + 1
    methylation_string = line[delim_loc[10]+1:delim_loc[11]] 
    mvec_string =  line[delim_loc[11]+1:delim_loc[12]] 
    bsseq_string =  line[delim_loc[12]+1: len(line) - 1] 
    #read_name_str = line[delim_loc[8] + 1:delim_loc[9]]
    #read_name = read_name_str[0:read_name_str.rindex("_")] 
    methylation_flank_string = methylation_string [m_vec_start:m_vec_stop]
    mvec_flank_string = mvec_string[m_vec_start:m_vec_stop]
    bsseq_flank_string = bsseq_string[m_vec_start:m_vec_stop]
    
    first_col = "`".join([line[delim_loc[2] + 1:delim_loc[3]], 
                          line[delim_loc[8] + 1:delim_loc[9]]])
                          #sam_flag_dict[read_name]]) 
    # record the real footprint lengths for footprints in the flank # 
    # example: `...........FFFFFFFFFFFFF.....FFF` 
    # ############## 
    to_print = "\t".join ([first_col, methylation_flank_string, str(m_vec_start), str(m_vec_stop), "strand" + strand, str(lflank), str(rflank)]) 
    #dbg_file.write(to_print + "\n")
    real_footprint_length_dict[first_col], footprint_length_at_bp_res_dict[first_col] = get_real_footprint_length(
         methylation_flank_string, m_vec_start, m_vec_stop, methylation_string)
    if strand == "+":
        to_write = "\t".join([first_col, methylation_flank_string, 
                              mvec_flank_string, bsseq_flank_string])
        output_fp.write(to_write + "\n") 
        footprint_dict[first_col] = methylation_flank_string 
        footprint_length_at_bp_res_tsv_fp.write(
           first_col + "\t" + "\t".join(map(str, 
           footprint_length_at_bp_res_dict[first_col])) + "\n")
    elif strand == "-": 
        to_write = "\t".join([first_col, 
                              methylation_flank_string[::-1],
                              mvec_flank_string[::-1],
                              bsseq_flank_string[::-1]]) 
        output_fp.write(to_write + "\n") 
        footprint_dict[first_col] = methylation_flank_string[::-1]
        footprint_length_at_bp_res_dict[first_col] = footprint_length_at_bp_res_dict[first_col][::-1]
        footprint_length_at_bp_res_tsv_fp.write(
           first_col + "\t" + "\t".join(map(str, 
           footprint_length_at_bp_res_dict[first_col][::-1])) + "\n")
        
		# since reads with both sam flags (99~147, and 83~163) are on 5'-3', I am
    # saving the footprint string "as is" also to plot both types of reads 
    # around the peak center
    to_write = "\t".join([first_col, methylation_flank_string])
    
    strand_agnostic_fp.write(to_write + "\n")
    strand_agnostic_footprint_dict[first_col] = methylation_flank_string 
    
pickle.dump(real_footprint_length_dict, real_length_dict)       
pickle.dump(footprint_dict, footprint_dict_fp)
pickle.dump(strand_agnostic_footprint_dict, strand_agnostic_footprint_dict_fp)
pickle.dump(footprint_length_at_bp_res_dict, footprint_length_at_bp_res_fp)

output_fp.close() 
input_fp.close()
real_length_dict.close()
footprint_dict_fp.close()
strand_agnostic_footprint_dict_fp.close() 
strand_agnostic_fp.close() 
footprint_length_at_bp_res_fp.close()
footprint_length_at_bp_res_tsv_fp.close()
#dbg_file.close()
