import os
import sys
import re
#SRR3133326.2_2/1	83	chr3R	3248440	40	150M	=	3248197	-393	ATACAAAAATATTATATAGAACTTTCTAAATCTCATAATTATATACAAAACCCACACAATTCTCAACAAAATCGAACAAAAAAAAAAAATTCTAACCAAAACCTTCAAACTGACCAAAAAAATAAAAAACCTTAACCTATACTTATTTAT	G.G.IIIGGAGGGG<GAGAGGAGGA<.GGA<...GAIGGGG.A.G.GGGGGG<..AGGGGGGAGGAAIIIGGGGGG.IIIIIIIIIIIIIIIIIIIIIIGGIIGGIIIGIIIIIGIIIIIIIIIIIIIGIIIIIGIIIIGGIIIIGGGGG	NM:i:38	MD:Z:2G3G1G1G6G2G6G1G6G0G6G3G0G5G1G6G0G3G7G0G0G4G0G1G5G0G2G2G5G1G3G4G5G2G6G4G1G6C2	XM:Z:..h...h.h.h......hH.h......x.h......hh......h...hh.....z.z......zx...h...Z...zxh....hh.h.....xh..z..h.....z.h..Xh....h.....h..h......h....x.h.........	XR:Z:CT	XG:Z:GA
#SRR3133326.3_3/1	99	chr2R	16094281	40	149M	=	16094425	294	AAGGGAAGATATTTATATTATTTTAAAAATGGTTAAAATATAGATGATATATTTTTTTTTTTGCTGTGTATAGGGTTATAAGTATATCGCTGTTGATGTTATTTATCGCAATTTTGCCGTTTGCTATTGCTGCAGTAGCAATATGGTCA	GGGGGGGGGIIIIGIIIIIIIIIIGGGGGGIGGIGIIIIGIIIGGGGIIIGGIIIIGIIIIIGGIGGGGIIGGGGIIGIIIGIGGGGGGGGGGIGGGGGGGGGGGGAGGIGGGIIAGGAGGI<GG<GG.<GGGGGGG.GGGIGGGGGII	NM:i:24	MD:Z:9C5C7C5C19C1C7C1C6C1C4C0C1C3C9C8C3C5C8C3C1C8C5C4C2	XM:Z:.........h.....h.......h.....z...................h.h.......h.z.X....h.x....hh.h...h....Z.X..x........h...xZ.H..h....XZ..x..Hh.x..X..X..x..H..h....hH.	XR:Z:CT	XG:Z:CT

############### Idea of suppressing all non CG/GC methylation ###########################
#
#
#    AAGGGAAGATATTTATATTATTTTAAAAATGGTTAAAATATAGATGATATATTTTTTTTTTTGCTGTGTATAGGGTTATAAGTATATCGCTGTTGATGTTATTTATCGCAATTTTGCCGTTTGCTATTGCTGCAGTAGCAATATGGTCA	
#    .........h.....h.......h.....z...................h.h.......h.z.X....h.x....hh.h...h....Z.X..x........h...xZ.H..h....XZ..x..Hh.x..X..X..x..H..h....hH.
#    
#    all Z, z, will be kept. All H, h, X, x will be checked for their upstream single base - if it is G it will be kept else will be replced with `.`.  
#
#########################################################################################


def suppress_non_cg_gc(bs_seq, mvec, pattern_type = None):
    idx_loc = []
    to_suppress = [] # filtered locations 
    mvec_len = len(mvec)
    suppressed_mvec = []
    for idx in range(mvec_len):
        if mvec[idx] in [".", "Z", "z"]:
            suppressed_mvec.append(mvec[idx])
        else:
            if pattern_type == "HCH": 
                if idx == 0:
                    suppressed_mvec.append(".")
                elif bs_seq[idx - 1] == "G":
                    suppressed_mvec.append(mvec[idx])
                else:
                    suppressed_mvec.append(".")
            elif pattern_type == "DGD":
                if idx == mvec_len - 1: 
                    suppressed_mvec.append(".")
                elif bs_seq[idx + 1] == "C":
                    suppressed_mvec.append(mvec[idx])
                else: 
                    suppressed_mvec.append(".")
    suppressed_mvec_str = "".join(suppressed_mvec)
    return suppressed_mvec_str 
          
        
for l in sys.stdin:
    d_loc = [m.start() for m in re.finditer("\t",l)] 
    sam_flag = l[d_loc[0] + 1 : d_loc[1]]
    bs_seq = l[d_loc[8] + 1: d_loc[9]]
    mvec = l[d_loc[12] + 1 + 5: d_loc[13]]
    if sam_flag in ["99", "147"]:
        suppressed_mvec_str = suppress_non_cg_gc(bs_seq, mvec, pattern_type = "HCH")
    elif sam_flag in ["83", "163"]:
        suppressed_mvec_str = suppress_non_cg_gc(bs_seq, mvec, pattern_type = "DGD")
    to_write = l[0:d_loc[12]] + "\t" + "XM:Z:" + suppressed_mvec_str + "\t" + l[d_loc[13] +1:len(l) - 1]
    print (to_write)
