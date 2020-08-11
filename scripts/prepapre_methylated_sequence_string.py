import os 
import sys
import re

mate_dict = {'83': '163', '99' : '147'}

###############################################################################
#                                Reference genome                             #
# 5' --------------------------------------------------------------------- 3' #
# 3' --------------------------------------------------------------------- 5' #
#                               |                                             #
#   Read                        |                                             #
# 5' ------------ 3' cond1      |                                             # 
#                  -----------> |                                             #
# 3' ------------ 5'            |                                             #
#             Mate              |                                             #
#                               |                                             #
#                               |                                             #
#                                 Mate --->        : sam flag 163             #
#                                 --------------                              #
#                                 --------------                              #
#                                        <--- Read : sam flag 83              #
#   Read                        |                                             #
# 5' ------------ 3' cond2      |                                             # 
#                  -----------> |                                             #
# 3' ------------ 5'            |                                             #
#             Mate              |                                             #
#                               |                                             #
#                               |                                             #
#                                 Read --->        : sam flag 99              #
#                                 --------------                              #
#                                 --------------                              #
#                                        <--- Mate : sam flag 147             #
###############################################################################

########################## why this code ######################################
# I am extracting sequences also in order to do random footprint analysis     #
#                                                                             #
#                                                                             #
#                                                                             #
#                                                                             #
###############################################################################
def get_fragment_methylation_vector (read1, read2):
    """
    Based on sam flag I am going to prepare methylation string. If two reads
    don't overlap then I am going to put `N`s in the gap.
    """
    read1_delim_loc = [m.start() for m in re.finditer("\t", read1)]
    read2_delim_loc = [m.start() for m in re.finditer("\t", read2)]
    read1_sam_flag = read1[read1_delim_loc[0]+1:read1_delim_loc[1]]
    read2_sam_flag = read2[read2_delim_loc[0]+1:read2_delim_loc[1]] 
    read1_start = int(read1[read1_delim_loc[2]+1:read1_delim_loc[3]])
    read2_start = int(read2[read2_delim_loc[2]+1:read2_delim_loc[3]])
    read1_len = read1_delim_loc[9] - read1_delim_loc[8] - 1 
    read2_len = read2_delim_loc[9] - read2_delim_loc[8] - 1
    complete_methylation_vec = None
    if read1_sam_flag in ['83', '99']: # thats the only possibilty we are considering
        read1_methylation_vector = read1[read1_delim_loc[8]+1: read1_delim_loc[9]]
        read2_methylation_vector = read2[read2_delim_loc[8]+1: read2_delim_loc[9]]
        chrom_name = read1[read1_delim_loc[1]+1:read1_delim_loc[2]]
        read_name = read1[0:read1_delim_loc[0]]
        
        if read1_sam_flag == '83': # read (5'-3'; first in pair) is on reverse strand 
            # meaning that the mate is present on forward strand and is in the left 
            # so assign read2 methylation vec in the start
            complete_methylation_vec = read2_methylation_vector
            if read2_start + read2_len - 1 >= read1_start: # read and mate overlap; 
                                                           # in this case 
                # oooooo
                # 012345
                #     oooooo
                #     456789 
                # read2_start + read2_len - read1_start = 0 + 6 - 4 = 2 
                start_pos_for_read1 =  read2_start + read2_len - read1_start 
                complete_methylation_vec = \
                  complete_methylation_vec + read1_methylation_vector[start_pos_for_read1:read1_len]
                bed_line = "\t".join([chrom_name, 
                                      str(read2_start), 
                                      str(read1_start + read1_len - 1),
                                      read_name + "_overlapping" + "`83~163",
                                      ".",
                                      complete_methylation_vec]) 
                print (bed_line)
            elif read2_start + read2_len - 1 < read1_start: 
                # oooooo
                # 012345
                #       oooo
                #       6789
                gap_len = read1_start - (read2_start + read2_len - 1) - 1
                number_Ns = 'N'*gap_len
                complete_methylation_vec = \
                   read2_methylation_vector + number_Ns + read1_methylation_vector
                bed_line = None
                if len(number_Ns) == 0: # adjacent
                    bed_line = "\t".join([chrom_name, 
                                          str(read2_start), 
                                          str(read1_start + read1_len - 1),
                                          read_name + "_adjacent" + "`83~163",
                                          ".",
                                          complete_methylation_vec]) 
                else: 
                    bed_line = "\t".join([chrom_name, 
                                          str(read2_start), 
                                          str(read1_start + read1_len - 1),
                                          read_name + "_gap" + "`83~163",
                                          ".",
                                          complete_methylation_vec]) 
                print (bed_line)
        elif read1_sam_flag == '99': # read (5'-3'; first in pair) is on forward strand 
            # meaning that the read is present on forward strand and is in the left 
            # so assign read1 methylation vec in the start
            complete_methylation_vec = read1_methylation_vector
            if read1_start + read1_len - 1 > read2_start: # read and mate overlap; 
                                                      # in this case 
                start_pos_for_read2 =  read1_start + read1_len - read2_start
                complete_methylation_vec = \
                  complete_methylation_vec + read2_methylation_vector[start_pos_for_read2:read2_len]
                bed_line = "\t".join([chrom_name, 
                                      str(read1_start), 
                                      str(read2_start + read2_len  - 1),
                                      read_name + "_overlapping" + "`99~147",
                                      ".",
                                      complete_methylation_vec]) 
                print (bed_line)
            elif read1_start + read1_len - 1 < read2_start:
                gap_len = read2_start - (read1_start + read1_len - 1) -  1
                number_Ns = 'N'*gap_len
                complete_methylation_vec = \
                   read1_methylation_vector + number_Ns + read2_methylation_vector
                bed_line = None
                if len(number_Ns) == 0: #adjacent
                    bed_line = "\t".join([chrom_name, 
                                          str(read1_start), 
                                          str(read2_start + read2_len - 1),
                                          read_name + "_adjacent" + "`99~147",
                                          ".",
                                          complete_methylation_vec]) 
                else:
                    bed_line = "\t".join([chrom_name, 
                                          str(read1_start), 
                                          str(read2_start + read2_len - 1),
                                          read_name + "_gap" + "`99~147",
                                          ".",
                                          complete_methylation_vec]) 
                print (bed_line)
                
    else:
        bed_line = "\t".join([chrom_name, 
                             str(read1_start), 
                             str(read2_start + read2_len - 1),
                             read_name + "_not_considered" + "`NA~NA",
                             ".",
                             "NA"]) 
        print (bed_line)
    
for l in sys.stdin:
    read1 = l
    read2 = next(sys.stdin)
    # split read line with `\t` as delimiter
  
    get_fragment_methylation_vector(read1, read2)
    # get the sam flag they have to correspond to ['83', '163', '99', '147']
