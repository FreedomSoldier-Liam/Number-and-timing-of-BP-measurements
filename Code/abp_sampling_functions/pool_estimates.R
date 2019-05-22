
# This is my rudimentary version
pool_estimates <- function(qhat, var_qhat, return.ci=FALSE){
  
  m = length(qhat)
  qbar = mean(qhat)
  ubar = mean(var_qhat)
  B = var(qhat)
  
  if(return.ci){
    list(
      mean = qbar,
      lower = qbar + qnorm(0.025) * sqrt(ubar + (1+1/m) * B),
      upper = qbar + qnorm(0.975) * sqrt(ubar + (1+1/m) * B)
    )
  } else {
    list(
      mean = qbar,
      variance = ubar + (1+1/m) * B
    )
    
  }
  
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# This version follow's Rubin's rules exactly
pool_estimates <- function(
  means,
  variances,
  nimputes,
  nobs,
  return.ci = FALSE
){
  vec = c()
  for (i in 1:nimputes)
    vec %<>% c(means[i], sqrt(variances[i]))
  
  out <- quiet(array(
    vec,
    dim = c(1, 2, nimputes),
    dimnames = list(1, c("Estimate", "Std.Error"), 1:nimputes)
  ) %>%
      miWQS::pool.mi(n = nobs))
  
  if (return.ci) {
    list(
      mean = as.numeric(out['pooled.mean']),
      variance = as.numeric(out['pooled.total.se']) ^ 2,
      lower = as.numeric(out['CI.1']),
      upper = as.numeric(out['CI.2'])
    )
  } else {
    list(mean = as.numeric(out['pooled.mean']),
      variance = as.numeric(out['pooled.total.se']) ^ 2)
  }
}

