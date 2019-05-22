
library(mice)

nimpute <- 10 
nboots <- 100

data_labels = c('JHS','CARDIA')

# winners object is created in part c
times = winners$overall %>% 
  stringr::str_replace_all(pattern='-', replacement = '_') %>% 
  purrr::set_names()


bootstrapper <- function(data_label, times, params){
  
  nobs <- file.path(
    "Datasets", paste(data_label,'analysis.RDS',sep='_')
  ) %>% 
    readRDS() %>% 
    nrow()
  
  map(times,
    .f = function(sampling_time){
      
      if(params$impute==TRUE){
        
        results <- seq(nimpute) %>% 
          purrr::set_names() %>% 
          map(
            ~ do_prevalence_comparison(
              data_label = data_label,
              times = sampling_time,
              nboots = nboots,
              params = params,
              abpm_impute_iter = .
            )
          )
        
      } else {
        
        results <- do_prevalence_comparison(
          data_label = data_label,
          times = sampling_time,
          nboots = nboots,
          params = params,
          abpm_impute_iter = 1
        )
        
      }
      
      difference_results <- results %>% 
        map(~.x$boot_diff) %>% 
        bind_rows(.id = 'impute') 
      
      estimate_results <- results %>% 
        map(~.x$boot_estm) %>% 
        bind_rows(.id = 'impute') 
      
      summary_results <- results %>% 
        map(~.x$boot_smry) %>% 
        bind_rows(.id = 'impute') %>% 
        pivot_wider(
          names_from = 'bp_type', 
          values_from = c('pr_mean', 'pr_var')
        ) %>% 
        group_by(outcome, bp_msr) %>% 
        nest() %>% 
        mutate(
          pooled_diffs = map(
            data,
            ~ pool_estimates(
              means = .$diff_mean,
              variances = .$diff_var,
              return.ci = TRUE,
              nobs = nobs, 
              nimputes = nimpute
            ) %>% 
              dplyr::bind_cols()%>%
              select(
                diff_mean = mean,
                diff_lwr = lower,
                diff_upr = upper
              )
          ),
          pooled_pr_tru = map(
            data,
            ~ pool_estimates(
              means = .$pr_mean_tru,
              variances = .$pr_var_tru,
              return.ci = TRUE,
              nobs = nobs, 
              nimputes = nimpute
            ) %>%
              dplyr::bind_cols() %>%
              select(
                pr_tru_mean = mean,
                pr_tru_lwr = lower,
                pr_tru_upr = upper
              )
          ),
          pooled_pr_est = map(
            data,
            ~ pool_estimates(
              means = .$pr_mean_est,
              variances = .$pr_var_est,
              return.ci = TRUE,
              nobs = nobs, 
              nimputes = nimpute
            ) %>% 
              dplyr::bind_cols() %>%
              select(
                pr_est_mean = mean,
                pr_est_lwr = lower,
                pr_est_upr = upper
              )
          )
        ) %>% 
        dplyr::select(outcome,bp_msr,starts_with("pooled")) %>% 
        unnest()
      
      list(
        diff = difference_results,
        estm = estimate_results,
        smry = summary_results
      )
      
    }
  )
}


# params object is created in markdown
prev_comp <- data_labels %>% 
  purrr::set_names() %>% 
  purrr::map( 
  bootstrapper, 
  times = times,
  params = params
)

saveRDS(prev_comp,file = 'Results/results_prev_comp.RDS')



