import os
import sys 
import re
# sys.argv[1]:  peak center bed followed by fragment bed
# sys.argv[2]: left flank 
# sys.argv[3]: right flank
# sys.argv[4]: output file

# example string to process: 

############ plus strand ###############
# chr2L   11806729        11806730        chr2L:11806729-11806729^20      .       +       chr2L   11806644        11806923        SRR3133326.203452_203452/1_overlapping  .       .h..h.H..x...h.H.Z.H.........Z....x...x........X...........h......h..zxhhhh.........z.....................xh...h.....hh..x...x..xhh.......h..h.h........h....h..z..........h.................X......Z.H.......h.....Zx....H.............h........Z....Z..........z.....X................
########################################

############ minus strand ###############
# chr2L   13205596        13205597        chr2L:13205596-13205596^20      .       -       chr2L   13205465        13205732        SRR3133326.98156_98156/1_overlapping    .       .hh.hh.h.............................hhh..............Z...hhh...Z...........x....x........hhx..Z......hh............h.h.h......xz.hx...z.....z..x..h.........h.....z....h.....h.............h................hh..z.xz....z......H...Z.H...hx.....h.Z..........Z...Z...Z..xZ.
########################################

input_fp = open(sys.argv[1])
lflank = int(sys.argv[2])
rflank = int(sys.argv[3])
out_fp = open(sys.argv[4], "w")
#dbg_fp = open ("tmp/first.tsv", "w")

for line in input_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    # check the strand first
    strand = line[d_loc[4] + 1: d_loc[5]]
    if strand == "+": 
        peak_start = int(line[d_loc[0] + 1: d_loc[1]])
        read_start = int(line[d_loc[6] + 1: d_loc[7]])
        read_end = int(line[d_loc[7] + 1: d_loc[8]])
        if (peak_start - read_start >= lflank) and \
           (read_end - peak_start >= rflank): 
           out_fp.write(line)
           #dbg_fp.write("\t".join([str(peak_start), str(read_start), str(read_end), str(peak_start - read_start), 
                                   #str(read_end - peak_start)] ) + "\n")
        else:
            continue
    elif strand == "-":
        peak_start = int(line[d_loc[0] + 1: d_loc[1]])
        read_start = int(line[d_loc[6] + 1: d_loc[7]])
        read_end = int(line[d_loc[7] + 1: d_loc[8]])
        if (peak_start - read_start >= rflank) and \
           (read_end - peak_start >= lflank): 
           out_fp.write(line)
           #dbg_fp.write("\t".join([str(peak_start), str(read_start), str(read_end), str(peak_start - read_start), 
                                   #str(read_end - peak_start)] ) + "\n")
        else:
            continue
    elif strand == ".": # closed enhancers peak strand is not available ; just copying the + strand code here
        peak_start = int(line[d_loc[0] + 1: d_loc[1]])
        read_start = int(line[d_loc[6] + 1: d_loc[7]])
        read_end = int(line[d_loc[7] + 1: d_loc[8]])
        if (peak_start - read_start >= lflank) and \
           (read_end - peak_start >= rflank): 
           out_fp.write(line)
           #dbg_fp.write("\t".join([str(peak_start), str(read_start), str(read_end), str(peak_start - read_start), 
                                   #str(read_end - peak_start)] ) + "\n")
        else:
            continue        
    else: 
        continue
out_fp.close()
input_fp.close()
#dbg_fp.close()
