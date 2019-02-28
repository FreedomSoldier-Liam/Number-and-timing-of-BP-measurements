gen_rowname.continuous_variable <- function(x){
  paste(x$label,x$unit,sep=', ')
}

gen_rowname.categorical_variable <- function(x){
  if(length(x$levels)>2){
    paste(x$label)
  } else {
    paste(x$levels[2])
  }
}