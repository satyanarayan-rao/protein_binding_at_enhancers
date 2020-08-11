import sys
import re
inp_fp = open(sys.argv[1])
gap = int(sys.argv[2])
wobble_string_list = ["F" + "".join (["."]*(i+1))  + "F" for i in range(gap) ]
out_fp = open(sys.argv[3], "w")

for line in inp_fp:
    # example of line:
    # 
    # chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147	.	......FFFFF...............FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFF..............FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFFFFFFFFFFFFFFFFFFFFFFFFFF.FFFFFFF.FFFFFFFFFFF....................................	..X..Z..x..Z.............H.........h...h.h..............................Hh..Z.......H....H........h....hh..h....hh....z....z...x....x....z....z.....z....h...h....H.............h.....h........x........z...H..h.hh.z..hx...x...h..h.x..H....h..H......h....Z......HH.h.hh.z...h.h.hhh.h.h.hh...	GGCTGCGTTTGCGTGGAAAGAAGAGCAAAAATATTTTAATTTTAATTGGAAAGGGAATGGGAATTGGAATTGCTTTCGGGGGGGCATTGCAAATTTTATATTTTTAATTATTTTAATTTGAAATGAGTAGTGTTGTGTGGTGTGGTTTTGGTGTTAATAAATCAAATTTATTGTATTTATATTTAGAAGAGTTTTATGTGTGTGCATTTTTTTGATTAGTTTGTTTTTTTTGCTTTTTTTCTAGATTTTAGACGGGGGGCCTTTTTTTGTTTTTTTTTTTTTATTATT
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    footprint_str = line[d_loc[10] + 1: d_loc[11]]
    footprint_str_list = list(footprint_str)
    for w in wobble_string_list:
        start = 0
        ret = 0
        while ret != -1:
            ret = footprint_str[start:].find(w)
            if ret != -1:
                for i in range(len(w)):
                    footprint_str_list[start + ret + i] = "F"
                start = start + ret + len(w) - 1
            elif ret == -1: 
                break
    to_write = line[0: d_loc[10]] + "\t" + "".join(footprint_str_list)  + "\t" + line[d_loc[11] + 1:]
            
    out_fp.write(to_write)

inp_fp.close()
out_fp.close()
