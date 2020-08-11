matrix_layout = function (reserved, count_vec){
    # 
    total = sum(count_vec)
    frac_vec = round (count_vec/total, 2)
    contib = (round ((1 - reserved)/sum(frac_vec), 2))*100
    dominance = frac_vec*contib 
    r_hundred = round(reserved,2)*100
    to_prep_matrix = c (r_hundred, dominance)
    z = c()
    cnt = 1
    for (i in to_prep_matrix){
        pseudo_i = i 
        if (pseudo_i < 3){
            pseudo_i = 3
        }
        z = c (z, rep(cnt, pseudo_i))
        cnt = cnt + 1 
    }
    lmatrix = matrix(z, length (z), 1, byrow = T)
    return(lmatrix)
}
