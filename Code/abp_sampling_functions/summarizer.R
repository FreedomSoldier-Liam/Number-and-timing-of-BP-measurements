

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
        ) 
      
      tlab=paste(sampling_times,collapse='_')
      tlab=gsub('.','deci',tlab,fixed = TRUE)
      fname=paste0(paste(input,collapse='_'),"_",tlab,".RDS")
      
      saveRDS(results,file.path('Results','Full Results',fname))
      
      if(!input$impute){
        
        results %<>% 
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
        
        sbp_mae <- results$slp_sbp_est %>% 
          subtract(results$slp_sbp) %>% 
          abs() %>% 
          boot_vec_mean()
        
        dbp_mae <- results$slp_dbp_est %>% 
          subtract(results$slp_dbp) %>% 
          abs()  %>% 
          boot_vec_mean()
        
        nht_classif <- (results$nht==results$nht_est)  %>% 
          boot_vec_mean()
        
        kappa <- Kappa.test(x=results$nht,y=results$nht_est)$Result
        
        kappa_variance <- kappa$conf.int %>% 
          diff() %>% 
          divide_by(2) %>% 
          divide_by(1.96) %>% 
          raise_to_power(2)
        
        output <- tibble(
          sbp_mae_est     = sbp_mae$mean,
          sbp_mae_var     = sbp_mae$variance,
          dbp_mae_est     = dbp_mae$mean,
          dbp_mae_var     = dbp_mae$variance,
          nht_classif_est = nht_classif$mean,
          nht_classif_var = nht_classif$variance,
          nht_kappa_est   = kappa$estimate,
          nht_kappa_var   = kappa_variance,
          label           = paste(sampling_times, collapse = '-')
        )
        
      } else {
        
        output_unpooled <- seq(1, 10) %>% 
          map(
            .f=function(iter){
              
              tmp = results %>%
                dplyr::mutate(
                  sbp = .[[paste0('sbp_imp_',iter)]],
                  dbp = .[[paste0('dbp_imp_',iter)]]
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
              
              sbp_mae <- tmp$slp_sbp_est %>% 
                subtract(tmp$slp_sbp) %>% 
                abs() %>% 
                boot_vec_mean()
              
              dbp_mae <- tmp$slp_dbp_est %>% 
                subtract(tmp$slp_dbp) %>% 
                abs()  %>% 
                boot_vec_mean()
              
              nht_classif <- (tmp$nht==tmp$nht_est)  %>% 
                boot_vec_mean()
              
              kappa <- Kappa.test(x=tmp$nht,y=tmp$nht_est)$Result
              
              kappa_variance <- kappa$conf.int %>% 
                diff() %>% 
                divide_by(2) %>% 
                divide_by(1.96) %>% 
                raise_to_power(2)
              
              tibble(
                sbp_mae_est     = sbp_mae$mean,
                sbp_mae_var     = sbp_mae$variance,
                dbp_mae_est     = dbp_mae$mean,
                dbp_mae_var     = dbp_mae$variance,
                nht_classif_est = nht_classif$mean,
                nht_classif_var = nht_classif$variance,
                nht_kappa_est   = kappa$estimate,
                nht_kappa_var   = kappa_variance,
                label           = paste(sampling_times, collapse = '-')
              )
            }
          ) %>% 
          bind_rows()
        
        sbp_pool <- pool_estimates(
          means = output_unpooled$sbp_mae_est,
          variances = output_unpooled$sbp_mae_var,
          nimputes = 10,
          nobs = nrow(data)
        )
        
        dbp_pool <- pool_estimates(
          means = output_unpooled$dbp_mae_est,
          variances = output_unpooled$dbp_mae_var,
          nimputes = 10,
          nobs = nrow(data)
        )
        
        nht_pool <- pool_estimates(
          means = output_unpooled$nht_classif_est,
          variances = output_unpooled$nht_classif_var,
          nimputes = 10,
          nobs = nrow(data)
        )
        
        kpa_pool <- pool_estimates(
          means = output_unpooled$nht_kappa_est,
          variances = output_unpooled$nht_kappa_var,
          nimputes = 10,
          nobs = nrow(data)
        )
        
        output <- tibble(
          sbp_mae_est     = sbp_pool$mean,
          sbp_mae_var     = sbp_pool$variance,
          dbp_mae_est     = dbp_pool$mean,
          dbp_mae_var     = dbp_pool$variance,
          nht_classif_est = nht_pool$mean,
          nht_classif_var = nht_pool$variance,
          nht_kappa_est   = kpa_pool$mean,
          nht_kappa_var   = kpa_pool$variance,
          label           = paste(sampling_times, collapse = '-')
        )
        
      }
      
      output
      
    }) %>% 
      bind_rows()
    
  }) %>% 
    bind_rows()
  
}
