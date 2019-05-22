
pr_comparison <- function(
  outcomes,
  exposures,
  control, 
  imputes
){
  
  names(outcomes)<- outcomes
  
  model_formulas <- outcomes %>% 
    map(
      .f = function(outcome){
        map(
          exposures, 
          .f = function(exposure){
            as.formula(
              paste(
                outcome, '~', paste(exposure,collapse='+'), "+",
                paste(
                  control, collapse = ' + '
                )
              )
            )
          }
        )
      }
    )
  
  model_estimates <- model_formulas %>% 
    # for each model formula, 
    map(
      .f = function(outcome_formulas){
        map(
          outcome_formulas,
          .f = function(model_formula){
            # fit poisson models to all multiply imputed datasets
            imputes %>% 
              map(
                .f = function(imputed_data){
                  imputed_data$id = 1:nrow(imputed_data)
                  geepack::geeglm(
                    formula = model_formula,
                    family = 'poisson',
                    id = id,
                    data = imputed_data
                  )
                }
              ) %>% 
              # pool the results
              mitml::testEstimates() %>%
              # pull out the estimates
              magrittr::use_series("estimates") %>%
              # coerce results into a tibble
              dplyr::as_tibble(rownames = "term")
          }
        ) %>% 
          # stack results together
          dplyr::bind_rows() %>% 
          dplyr::filter(term %in% flatten_chr(exposures)) %>% 
          dplyr::select(
            term,
            estimate = Estimate,
            std.error = Std.Error,
            p.value = `P(>|t|)`
          )
      }
    ) %>% 
    dplyr::bind_rows(.id='outcome')%>% 
    dplyr::group_by(outcome) %>%
    dplyr::mutate(pr = exp(estimate)) %>%
    separate('term', into = c('junk','bp_msr','bp_type')) %>% 
    select(-junk)
  
  pr_difference = model_estimates %>%
    arrange(outcome, bp_msr) %>% 
    group_by(outcome, bp_msr) %>% 
    dplyr::summarise(
      pr_diff = diff(pr),
      est_diff = diff(estimate)
    )
  
  list(
    model_estimates = model_estimates,
    pr_difference = pr_difference
  )
  
}

do_prevalence_comparison <- function(
  data_label,
  times,
  nboots,
  params,
  abpm_impute_iter
){
  
  times=gsub(".","deci",times, fixed=TRUE)
  
  control = if(data_label == 'CARDIA'){
    c('age',
      'sex',
      'race',
      'edu',
      'diabetes',
      'currentsmoker',
      'bpmeds',
      'slp_duration'
      )
  } else {
    c('age',
      'sex',
      'edu',
      'diabetes',
      'currentsmoker',
      'bpmeds',
      'slp_duration'
      )
  }
  
  analysis <- file.path(
    "Datasets", paste(data_label,'analysis.RDS',sep='_')
  ) %>% 
    readRDS()
  
  if(params$impute==TRUE){
    abpm_impute_names <- paste(
      c("sbp_imp","dbp_imp"),abpm_impute_iter,sep="_"
    )
  } else {
    abpm_impute_names <- paste(
      c("sbp","dbp")
    )
  }
  
  # Data import -------------------------------------------------------------
  
  full_result <- file.path(
    "Results",
    "Full Results",
    paste0(
      data_label,"_",
      params$impute,"_",
      params$strata,"_",
      params$time_variable,"_",
      paste(times, collapse='_'),
      '.RDS'
    )
  ) %>% 
    read_rds() %>% 
    group_by(subjid) %>% 
    dplyr::select(
      subjid,
      sbp_imp = abpm_impute_names[1],
      dbp_imp = abpm_impute_names[2]
    ) %>% 
    dplyr::summarize(
      slp_sbp_est=mean(sbp_imp),
      slp_dbp_est=mean(dbp_imp)
    ) %>% 
    mutate(
      nht_est = ifelse(
        slp_sbp_est >= 120 | slp_dbp_est >= 70, 1, 0
      )
    ) %>% 
    left_join(analysis, by = 'subjid') %>% 
    dplyr::rename(
      slp_sbp_tru = slp_sbp,
      slp_dbp_tru = slp_dbp
    ) %>% 
    mutate(
      slp_sbp_tru = slp_sbp_tru/10,
      slp_sbp_est = slp_sbp_est/10,
      slp_dbp_tru = slp_dbp_tru/10,
      slp_dbp_est = slp_dbp_est/10
    )
  
  # Bootstrapping ---------------------------------------------------------
  
  exposures <- list(
    true = c('slp_sbp_tru','slp_dbp_tru'),
    estm = c('slp_sbp_est','slp_dbp_est')
  )
  
  outcomes <- c("lvh","albuminuria")
  
  boot_results <- seq(nboots) %>% 
    set_names(paste(.)) %>% 
    map(
      .f=function(boot_replicate){
        
        set.seed(boot_replicate)
        
        boot_data <- full_result %>% 
          .[sample(1:nrow(.), replace=TRUE), ] %>% 
          droplevels()
        
        imputes <- boot_data %>% 
          mutate_if(is.character, as.factor) %>% 
          dplyr::select(-subjid) %>% 
          mice::mice(method='pmm',m = nimpute) %>% 
          mice::complete(action = 'all') %>% 
          map(as_tibble) %>% 
          map(mutate_at, outcomes, ~as.numeric(.x)-1)
        
        pr_comparison(
          outcomes = outcomes,
          exposures = exposures,
          control = control,
          imputes = imputes
        )
      }
    ) 
  
  boot_estimates <- boot_results %>% 
    map(~.$model_estimates) %>% 
    dplyr::bind_rows(.id='replication')
  
  boot_pr_diffs <- boot_results %>% 
    map(~.$pr_difference) %>% 
    dplyr::bind_rows(.id='replication')
  
  estimates <- boot_estimates %>% 
    dplyr::group_by(outcome, bp_msr, bp_type) %>% 
    dplyr::summarise(
      pr_mean = mean(pr),
      pr_var = var(pr)
    )
  
  pr_diffs <- boot_pr_diffs %>% 
    group_by(outcome, bp_msr) %>% 
    summarise(
      diff_mean = mean(pr_diff),
      diff_var = var(pr_diff)
    )
  
  list(
    boot_estm = boot_estimates,
    boot_diff = boot_pr_diffs,
    boot_smry = left_join(estimates, pr_diffs)
  )
  
}
