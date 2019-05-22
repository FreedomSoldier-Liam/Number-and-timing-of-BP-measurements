
library(kableExtra)
library(tidyverse)
library(magrittr)
library(gt)

files <- list.files(
  path      = 'Results',
  all.files = TRUE,
  full.names= TRUE, 
  pattern   = paste(
    params$impute,
    params$strata,
    params$time_variable,
    sep='_'
  )
)

jhs_indx <- map_lgl(files, grepl, pattern='JHS', fixed=T)
cardia_indx <- map_lgl(files, grepl, pattern='CARDIA', fixed=T)

names(files)[jhs_indx] <- 'JHS'
names(files)[cardia_indx] <- 'CARDIA'

results <- map_dfr(files, readRDS, .id='dataset') %>% 
  mutate(nmeasures = 1 + str_count(label, pattern='-')) %>% 
  arrange(dataset, protocol, nmeasures) 

winners <- list(
  overall = results %>% 
    group_by(protocol, nmeasures, label) %>% 
    summarise(mean_kappa = mean(nht_kappa_est)) %>% 
    arrange(protocol, nmeasures, desc(mean_kappa)) %>% 
    group_by(protocol, nmeasures) %>% 
    mutate(rank = 1:n()) %>% 
    filter(rank == 1) %>% 
    pluck("label"),
  bystudy = results %>% 
    group_by(protocol, nmeasures, dataset) %>% 
    arrange(desc(nht_kappa_est)) %>%
    slice(1) %>% 
    group_by(protocol, nmeasures) %>% 
    pluck("label") 
)

generate_sample_comparison_table <- function(
  winners, data, params
){
  data %>% filter(label %in% winners) %>% 
    mutate(
      sbp_mae = paste_ci(
        sbp_mae_est,
        sbp_mae_lwr,
        sbp_mae_upr
      ),
      dbp_mae = paste_ci(
        dbp_mae_est,
        dbp_mae_lwr,
        dbp_mae_upr
      ),
      nht_classif = paste_ci(
        nht_classif_est,
        nht_classif_lwr,
        nht_classif_upr
      ),
      nht_kappa = paste_ci(
        nht_kappa_est,
        nht_kappa_lwr,
        nht_kappa_upr
      ),
      label = gsub("-", ", ", label)
      #,label = gsub(".5", ":30", label)
    ) %>% 
    dplyr::select(
      dataset, protocol,
      `Sampling protocol` = label,
      `Mean absolute error (95% CI)_Asleep Systolic BP` = sbp_mae,
      `Mean absolute error (95% CI)_Asleep Diastolic BP` = dbp_mae,
      `Classification accuracy (95% CI)_% correctly classified` = nht_classif,
      `Classification accuracy (95% CI)_Kappa statistic` = nht_kappa
    ) %>% 
    dplyr::arrange(protocol,dataset) %>% 
    group_by(protocol,dataset) %>% 
    gt::gt(
      stub_group.sep = ", "
    ) %>% 
    gt::tab_header(
      title = 'Mean absolute error and classification accuracy of different sampling protocols.'
    ) %>% 
    gt::cols_split_delim(
      delim = '_',
      columns = c(
        'Mean absolute error (95% CI)_Asleep Systolic BP',
        'Mean absolute error (95% CI)_Asleep Diastolic BP'
      )
    ) %>% 
    gt::cols_split_delim(
      delim = '_',
      columns = c(
        'Classification accuracy (95% CI)_% correctly classified',
        'Classification accuracy (95% CI)_Kappa statistic'
      )
    ) %>% 
    gt::cols_align(
      columns = c(
        'Sampling protocol'
      ),
      align = 'center'
    ) %>% 
    gt::cols_align(
      columns = c(
        'Mean absolute error (95% CI)_Asleep Systolic BP',
        'Mean absolute error (95% CI)_Asleep Diastolic BP',
        'Classification accuracy (95% CI)_% correctly classified',
        'Classification accuracy (95% CI)_Kappa statistic'
      ),
      align = 'center'
    ) %>% 
    gt::tab_footnote(
      footnote = c(
        "Distributed protocol: A number of blood pressure measurements are taken between 1am and 5am, with at least 1 hour between each measurement",
        "Concentrated protocol: A number of blood pressure measurements are taken sequentially, with at most 30 minutes between each measurement, starting somewhere between 1am and 4am"
      ),
      locations = cells_column_labels(
        columns = c(
          'Sampling protocol'
        )
      )
    ) %>% 
    gt::tab_footnote(
      footnote = ifelse(
        params$impute==TRUE,
        "Ambulatory blood pressure measurements were adjusted by a suitable error term to account for oversampling.",
        "Ambulatory blood pressure measurements were not adjusted prior to sampling"
      ),
      locations = cells_title(groups='title')
    ) %>% 
    gt::tab_footnote(
      footnote = ifelse(
        params$time_variable=='tss',
        "blood pressure measurements were taken relative to the time that participants fell asleep",
        "blood pressure measurements were taken relative to time since midnight"
      ),
      locations = cells_title(groups='title')
    )
}



