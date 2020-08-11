from collections import defaultdict
def get_potential_cytosine_dict(): 
    """
    Returns a True/False defaultdict of cytosine context that could be methylated
    """
    potential_cytosines = defaultdict (lambda: defaultdict(lambda : False)) 
    sam_flag_99_147_list = ["CG", # {z|Z} CpG context 
               "CAG", "CCG", "CTG", # {x|X} CHG context 
               "CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT" # {h|H} CHH context
               ] # in `99~147` case we have to find these patterns and record the first position (as C is on the forward strand - so all starting Cs are targets) 
    sam_flag_83_163_list = ["CG", # {z|Z} CpG context 
               "CTG", "CGG", "CAG", # {x|X} CHG context 
               "TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG" # {h|H} CHH context
               ] # in `83~163` case we have to find these patterns and record the last position (as C is on the reverse strand - so all ending Gs are targets)
    for s in sam_flag_99_147_list: 
        potential_cytosines["99~147"][s] = True 
    for s in sam_flag_83_163_list: 
        potential_cytosines["83~163"][s] = True     
    return potential_cytosines 
