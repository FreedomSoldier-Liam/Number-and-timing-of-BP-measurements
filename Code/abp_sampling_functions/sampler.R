

abp_sampler<-function(data, times, time_variable, quality_check=FALSE){
  
  indx <- map_int(
    times, .f=function(sample_time){
      which.min(
        abs(data[[time_variable]]-sample_time)
      )
    } 
  )
  
  bad_sample <- any(data[[time_variable]][indx] - times > 1/2)
  
  if(quality_check & bad_sample){
    return(FALSE)
  } else if(quality_check & !bad_sample){
    return(TRUE)
  } else {
    if(bad_sample){
      return(NULL)
    } else {
      return(data[indx,])
    }
  }
}


boot_vec_mean <- function(x, xname, nboots = 1000, return.ci=FALSE){
  
  boot_results <- seq(1, nboots) %>% 
    map_dbl(
      ~ mean(x[sample(1:length(x),replace=TRUE)])
    )
  
  if(return.ci){
    list(
      median = median(boot_results),
      mean = mean(boot_results),
      lower = quantile(boot_results, 0.025, names = FALSE),
      upper = quantile(boot_results, 0.975, names = FALSE)
    )
  } else {
    
    list(
      mean = mean(boot_results),
      variance = var(boot_results)
    )
    
  }
  
  
}


