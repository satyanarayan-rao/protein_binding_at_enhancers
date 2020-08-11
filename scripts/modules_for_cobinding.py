def count_m_at_start_end(vec):
    len_vec = len(vec)
    st_cnt = 0 
    en_cnt = 0
    for i in range(len_vec):
        if vec[i] == "M":
            st_cnt +=1 
        else:
            break
    for i in range(len_vec - 1, -1, -1):
        if vec[i] == "M":
            en_cnt +=1 
        else:
            break
    return st_cnt, en_cnt 
def add_edge_colors(vec, st, en):
    vec_l = list(vec)
    for i in range (st, en + 1):
        vec_l[i] = 'E'
    ret_vec = "".join(vec_l)
    return ret_vec
def get_per_orange_for_each_footprint (fvec):
    flen_list = [] 
    cnt = 0 
    f_switch = False
    fvec_len = len(fvec)
    for c in fvec:
        if c == 'F':
            cnt +=1
            f_switch = True
        else:
            if cnt!=0:
            	flen_list.append(round((cnt/fvec_len)*100, 2))
            	cnt = 0 
    if cnt!= 0:
        flen_list.append(round((cnt/fvec_len)*100, 2))
    return flen_list
def get_read_start_and_end_from_lex_rex_reads(fvec):
    """
     0 1 2 3 4 5 6 7 8 9101112  
     M M F . . . F F F F M M M  
     return: 2, 9 
    """
    read_start = None
    read_end = None
    for i in range(len(fvec)):
        if fvec[i] != 'M':
            read_start = i
            break
    for i in range(len(fvec) - 1, -1, -1):
        if fvec[i] !='M':
            read_end = i
            break
    return read_start, read_end

def get_count_and_percentage_methylation (m_vec):
    total = 0
    total_lower = 0
    total_upper = 0
    per_c = 0
    for c in m_vec:
        if c == ".":
            continue
        elif c == c.lower(): 
            total +=1
            total_lower +=1
        elif c == c.upper():
            total +=1 
            total_upper +=1 
    if total>0:
        per_c = str(round(total_upper/total, 3)*100)
    else:
        per_c = "NA"
    return per_c, total_upper, total_lower, total

