import sys
inp_fp = open(sys.argv[1])
header = inp_fp.readline()

for line in inp_fp:
    l_items = line.strip().split()
    name = l_items[1]
    chrom = l_items[1].split(":")[0]
    st = l_items[1].split(":")[1].split("-")[0]
    en = l_items[1].split(":")[1].split("-")[1].split("^")[0]
    cl = l_items[1].split("^")[-1]
    strand = l_items[-1]
    to_write = "\t".join([chrom, st, en, name, ".", strand])
    print(to_write)
    
