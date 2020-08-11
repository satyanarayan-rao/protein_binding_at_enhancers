import os
import sys
import pickle

# sys.argv[1] = peak intersected with bam file
# sys.argv[2] = methylation status dictionary
# sys.argv[3] = output methylation matrix 
# sys.argv[4] = slop from peak center
input_bed = open(sys.argv[1])
methylation_status_dict = pickle.load(open(sys.argv[2], "rb")) 
output_fp = open (sys.argv[3], "w")
slop = int(sys.argv[4])
cnt = 0 
for line in input_bed: 
    line_items = line.strip().split() 
    chr_name = line_items[0]
    read_id = line_items[3]
    read_start = int(line_items[1])
    read_end = int (line_items[2])
    peak_start = int(line_items[13])
    peak_end = int(line_items[14])
    peak_center = int( (peak_start + peak_end)/2)
    flank_start = max (read_start, peak_center - slop) 
    flank_end = peak_center + slop
    # I don't care about flank end as it is always going to be positive 
    methylation_status_list = [line_items[15]]
    for chr_loc in range(flank_start, flank_end + 1): 
        dict_key = ":".join([read_id, chr_name, str(chr_loc)])
        if cnt < 5:
            print(dict_key)
        cnt += 1
        if dict_key in methylation_status_dict: 
            methylation_status_list.append(methylation_status_dict[dict_key])
        else:
            methylation_status_list.append("NA")

    to_write = "\t".join(methylation_status_list)
    output_fp.write(to_write + "\n")

output_fp.close() 
input_bed.close()   
