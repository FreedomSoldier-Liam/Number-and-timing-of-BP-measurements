
list_elements = function(x, 
                         end_string = 'and ', 
                         sep_string = '=',
                         col_string = ', ',
                         use_names=F){
  
  if(is.list(x)) x %<>% reduce(c)
  
  l = length(x)
  if(use_names & is.numeric(x)){
    x = paste(names(x), format(round(x,2),nsmall=2), sep = sep_string)
  } else if (use_names & is.character(x)){
    x = paste(names(x), x, sep = sep_string)
  }
  
  if (l == 1){
    paste(x)
  } else if (l == 2){
    paste0(x[1], ' and ', x[2])
  } else if (l >= 3){
    start = paste(x[1:(l-1)], collapse = col_string)
    stop =  paste0(end_string, x[l])
    paste(start, stop)
  }
}