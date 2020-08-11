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


def get_fragment_methylation_vector (read1, read2):
    """
    Based on sam flag I am going to prepare methylation string. If two reads
    don't overlap then I am going to put `N`s in the gap. 
    I am going to return bsseq string also. 
    Example string is shown below
    """
# SRR3133326.2_2/1	83	chr3R	3248440	40	150M	=	3248197	-393	ATACAAAAATATTATATAGAACTTTCTAAATCTCATAATTATATACAAAACCCACACAATTCTCAACAAAATCGAACAAAAAAAAAAAATTCTAACCAAAACCTTCAAACTGACCAAAAAAATAAAAAACCTTAACCTATACTTATTTAT	G.G.IIIGGAGGGG<GAGAGGAGGA<.GGA<...GAIGGGG.A.G.GGGGGG<..AGGGGGGAGGAAIIIGGGGGG.IIIIIIIIIIIIIIIIIIIIIIGGIIGGIIIGIIIIIGIIIIIIIIIIIIIGIIIIIGIIIIGGIIIIGGGGG	NM:i:38	MD:Z:2G3G1G1G6G2G6G1G6G0G6G3G0G5G1G6G0G3G7G0G0G4G0G1G5G0G2G2G5G1G3G4G5G2G6G4G1G6C2	XM:Z:..h...h.h.h......hH.h......x.h......hh......h...hh.....z.z......zx...h...Z...zxh....hh.h.....xh..z..h.....z.h..Xh....h.....h..h......h....x.h.........	XR:Z:CT	XG:Z:GA
# SRR3133326.2_2/1	163	chr3R	3248197	40	21M1I128M	=	3248440	393	TATAATATCCTTAATCCAAAATTTTTTTTTTATTTTAAAAAACAAATCATCCAACTAAAAAATCATATACAATTTACCAAATATCACATTAAACAAACCTATCTACTACACAAATTAACCCAAAAACATCTAAATTATAACTCATTTAAT	AGGGGIGGGGGGIIIIIIIIIIIGIIIIGIIGGGGIGIIGIIIGIIIIIIIIIIIGIIGGGGGGIIIGGGIGIIGIIGIIGGAGGGGIGG.AGGIGIIIIIGGGGGGIIGGGGGIGGGGGGAGGGGGGGGGIIGIAGGGIGIIIGGGGGG	NM:i:31	MD:Z:1G10G7G15G0G6G2G3G3G0G6G3G1G0G3G2G1G9G3G0G0G10G3G4G0G5G1G5G4G10G2	XM:Z:.h..........h.......h................hh......h..z...z...xh......z...h.zx...h..z.h.........h...zxh..........x...z....hh.....h.h.....x....h..........h..	XR:Z:GA	XG:Z:GA

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
        read1_methylation_vector = read1[read1_delim_loc[12]+6: read1_delim_loc[13]]
        read2_methylation_vector = read2[read2_delim_loc[12]+6: read2_delim_loc[13]] 
        read1_bsseq_vector = read1[read1_delim_loc[8]+1:read1_delim_loc[9]]
        read2_bsseq_vector = read2[read2_delim_loc[8]+1:read2_delim_loc[9]] 

        chrom_name = read1[read1_delim_loc[1]+1:read1_delim_loc[2]]
        read_name = read1[0:read1_delim_loc[0]]
        
        if read1_sam_flag == '83': # read (5'-3'; first in pair) is on reverse strand 
            # meaning that the mate is present on forward strand and is in the left 
            # so assign read2 methylation vec in the start
            complete_methylation_vec = read2_methylation_vector
            complete_bsseq_vec = read2_bsseq_vector
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
                complete_bsseq_vec = \
                  complete_bsseq_vec + read1_bsseq_vector[start_pos_for_read1:read1_len]
                bed_line = "\t".join([chrom_name, 
                                      str(read2_start), 
                                      str(read1_start + read1_len - 1),
                                      read_name + "_overlapping" + "`83~163",
                                      ".",
                                      complete_methylation_vec,
                                      complete_bsseq_vec]) 
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
                complete_bsseq_vec = \
                   read2_bsseq_vector + number_Ns + read1_bsseq_vector
                bed_line = None
                if len(number_Ns) == 0: # adjacent
                    bed_line = "\t".join([chrom_name, 
                                          str(read2_start), 
                                          str(read1_start + read1_len - 1),
                                          read_name + "_adjacent" + "`83~163",
                                          ".",
                                          complete_methylation_vec,
                                          complete_bsseq_vec]) 
                else: 
                    bed_line = "\t".join([chrom_name, 
                                          str(read2_start), 
                                          str(read1_start + read1_len - 1),
                                          read_name + "_gap" + "`83~163",
                                          ".",
                                          complete_methylation_vec,
                                          complete_bsseq_vec]) 
                print (bed_line)
        elif read1_sam_flag == '99': # read (5'-3'; first in pair) is on forward strand 
            # meaning that the read is present on forward strand and is in the left 
            # so assign read1 methylation vec in the start
            complete_methylation_vec = read1_methylation_vector
            complete_bsseq_vec = read1_bsseq_vector
            if read1_start + read1_len - 1 > read2_start: # read and mate overlap; 
                                                      # in this case 
                start_pos_for_read2 =  read1_start + read1_len - read2_start
                complete_methylation_vec = \
                  complete_methylation_vec + read2_methylation_vector[start_pos_for_read2:read2_len]
                complete_bsseq_vec = \
                  complete_bsseq_vec + read2_bsseq_vector[start_pos_for_read2:read2_len]
                bed_line = "\t".join([chrom_name, 
                                      str(read1_start), 
                                      str(read2_start + read2_len  - 1),
                                      read_name + "_overlapping" + "`99~147",
                                      ".",
                                      complete_methylation_vec, 
                                      complete_bsseq_vec]) 
                print (bed_line)
            elif read1_start + read1_len - 1 < read2_start:
                gap_len = read2_start - (read1_start + read1_len - 1) -  1
                number_Ns = 'N'*gap_len
                complete_methylation_vec = \
                   read1_methylation_vector + number_Ns + read2_methylation_vector
                complete_bsseq_vec = \
                   read1_bsseq_vector + number_Ns + read2_bsseq_vector
                bed_line = None
                if len(number_Ns) == 0: #adjacent
                    bed_line = "\t".join([chrom_name, 
                                          str(read1_start), 
                                          str(read2_start + read2_len - 1),
                                          read_name + "_adjacent" + "`99~147",
                                          ".",
                                          complete_methylation_vec,
                                          complete_bsseq_vec]) 
                else:
                    bed_line = "\t".join([chrom_name, 
                                          str(read1_start), 
                                          str(read2_start + read2_len - 1),
                                          read_name + "_gap" + "`99~147",
                                          ".",
                                          complete_methylation_vec,
                                          complete_bsseq_vec]) 
                print (bed_line)
                
    else:
        bed_line = "\t".join([chrom_name, 
                             str(read1_start), 
                             str(read2_start + read2_len - 1),
                             read_name + "_not_considered" + "`NA~NA",
                             ".",
                             "NA",
                             "NA"]) 
        print (bed_line)
    
for l in sys.stdin:
    read1 = l
    read2 = next(sys.stdin)
    # split read line with `\t` as delimiter
  
    get_fragment_methylation_vector(read1, read2)
    # get the sam flag they have to correspond to ['83', '163', '99', '147']
