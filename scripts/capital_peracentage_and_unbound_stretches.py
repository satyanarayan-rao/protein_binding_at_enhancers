def capital_percentage_and_stretch_of_unbound(mvec):
    len_vec = len(mvec)
    check_vec = "."*len_vec
    total_capital = 0 
    total_small = 0
    left_extreme = 0 
    right_extreme = len_vec 
    first_found = False
    base_counter = -1
    last_capital = -1   
    for c in mvec: 
        base_counter +=1
        if c == ".":
            continue 
        elif c == c.lower(): 
            total_small +=1 
        elif c == c.upper():
             total_capital +=1 
             if first_found == False: 
                 left_extreme = base_counter 
                 first_found = True
             else:
                 last_capital = base_counter
        else:
            continue 
    naked_dna_stretches = [] 
    if first_found == True: # there was at least one methylation event
         if last_capital == -1: # there was only one capital
             naked_dna_stretches.append(left_extreme) 
             naked_dna_stretches.append(0)
             naked_dna_stretches.append (right_extreme - left_extreme)
         else:
             naked_dna_stretches.append(left_extreme) 
             naked_dna_stretches.append(last_capital - left_extreme - 1)
             naked_dna_stretches.append(right_extreme - last_capital - 1)  
    else:
        naked_dna_stretches = [0, 0, 0]
    total = total_capital + total_small 
    percentage_capital = "NA"
    if total > 0:
    	percentage_capital = round(total_capital/total, 2)
    
    return len_vec, total, percentage_capital, naked_dna_stretches 
