

summarizer <- function(sampling_protocols, data, input, truth)
{
  map(sampling_protocols,.f=function(sampling_protocol){
    
    map(sampling_protocol, .f=function(sampling_times){
      
      results = data %>%
        plyr::ddply(
          .variables = 'subjid',
          .fun = abp_sampler, 
          times = sampling_times,
          time_variable = input$tvarbl
        ) %>% 
        group_by(subjid) %>% 
        dplyr::summarise(
          slp_sbp_est = mean(sbp),
          slp_dbp_est = mean(dbp)
        ) %>%
        mutate(
          nht_est=ifelse(
            slp_sbp_est >= 120 | slp_dbp_est >=70 , 1, 0
          )
        ) %>% 
        left_join(truth, by='subjid')
      
      tlab=paste(sampling_times,collapse='_')
      tlab=gsub('.','deci',tlab,fixed = TRUE)
      
      fname=paste0(paste(input,collapse='_'),"_",tlab,".RDS")
      
      saveRDS(results,file.path('Results','Full Results',fname))
      
      tibble(
        sbp_mae=mae(results$slp_sbp_est,results$slp_sbp),
        dbp_mae=mae(results$slp_dbp_est,results$slp_dbp),
        #lvm=cor(results$slp_sbp_est,results$lvm),
        nht_classif=mean(results$nht==results$nht_est),
        nht_kappa=Kappa.test(x=results$nht,y=results$nht_est)$Result$estimate,
        label=paste(sampling_times,collapse='-')
      )
      
    }) %>% 
      bind_rows()
    
  }) %>% 
    bind_rows()
  
}
