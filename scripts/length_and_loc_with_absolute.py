def get_real_footprint_length_with_abs_start (m_vec, m_vec_start, m_vec_stop, complete_vec):
    """
    extract complete footprint length from the selected methylation vector
    exmaple:
    complete vec: . . . . F F F F . . . F F F . . . . . F F F F F
    index       : 0 1 2 3 4 5 6 7 8 9 1011121314151617181920212223
    flank vec:  :             F F . . . F F F . . . . . F
                              6 7 8 9 10111213141516171819
    output: [4, 3, 5]
    """
    # first get all footprint lengths in m_vec
    flen_lengths = [len(a) for a in m_vec.split(".")]
    abs_loc = []
    if m_vec.startswith("F"):
        flen_lengths[0] = flen_lengths[0] +  len(complete_vec[0:m_vec_start].split(".")[-1])
        start_loc_first = m_vec_start
        for j in range(m_vec_start, -1, -1):
            if complete_vec[j] == 'F':
                start_loc_first -=1
            elif complete_vec[j] == '.':
                break
        abs_loc.append(start_loc_first+1)
    if m_vec.endswith("F") and len(flen_lengths)>0:
        flen_lengths[-1] = flen_lengths[-1] +  len(complete_vec[m_vec_stop:].split(".")[0])


    # exclude all zeros
    return_list = [a for a in flen_lengths if a!=0]
    # prepare a footprint length vector - at each index it will tell what is the footprint size that index it associated with
    first_f = False
    cnt = 0
    loc_first = []
    gap = True
    for c in m_vec:
        if c == 'F':
            first_f = True
            if gap == True:
                loc_first.append (cnt)
                gap = False
        else:
            first_f = False
            gap = True
        cnt +=1
    out_vec = [0]*len(m_vec)
    for i in range (len(loc_first)):
        for j in range(loc_first[i], len(m_vec)):
            if m_vec[j] == "F":
                out_vec[j] = return_list[i]
            else:
                break
    if len(loc_first) > 0:
        for i in loc_first:
            if i != 0: 
                abs_loc.append(m_vec_start + i)
    return return_list, out_vec, loc_first, abs_loc
