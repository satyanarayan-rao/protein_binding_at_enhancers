convert_string_to_numeric_vec = function(s){
  vec = c()
  s_c = unlist(strsplit(s, split = ""))
  for (c in s_c){
    if (identical(c, ".")){
      vec = c (vec, NA)
    }else if (identical(c, tolower(c))){
      vec = c (vec, 0)
    }else if (identical(c, toupper(c))){
      vec = c(vec, 1)
    }
  }
  return (vec)
}


convert_string_to_numeric_vec_zero_for_dot = function(s){
  vec = c()
  s_c = unlist(strsplit(s, split = ""))
  for (c in s_c){
    if (identical(c, ".")){
      vec = c (vec, 0)
    }else if (identical(c, tolower(c))){
      vec = c (vec, 0)
    }else if (identical(c, toupper(c))){
      vec = c(vec, 1)
    }
  }
  #print (vec)
  return (vec)
}
