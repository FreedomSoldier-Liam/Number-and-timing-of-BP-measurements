
paste_ci <- function(est, lwr, upr){
  
  paste0(
    adapt_round(est), ' (',
    adapt_round(lwr), ', ',
    adapt_round(upr), ')'
  )
  
}
