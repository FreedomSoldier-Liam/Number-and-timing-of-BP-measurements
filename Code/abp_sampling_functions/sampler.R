

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