import os
import sys
import re
import pickle 
from collections import defaultdict
from collections import OrderedDict
from extend_a_given_read import extend_footprint_from_center

footprint_fp = open(sys.argv[1]) 
footprint_dict = pickle.load(open(sys.argv[2], "rb")) 

for line in footprint_fp: 
    ex_footp, ex_mvec, ex_bsseq = extend_footprint_from_center(line, 
              footprint_dict, 0, 0, 15, 30, "+") 
    print ([ex_footp, ex_mvec, ex_bsseq])
