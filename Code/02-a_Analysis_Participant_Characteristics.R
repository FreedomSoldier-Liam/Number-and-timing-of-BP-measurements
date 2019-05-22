
suppressPackageStartupMessages(
  {
    library(kableExtra)
    library(tidyverse)
    library(magrittr)
    library(KableOne)
    library(labelled)
  }
)

read_rds(path = 'Datasets/pooled_analysis.RDS') %>% 
  select(-subjid) %>%
  tibble_one(strat = 'study', include.pval = TRUE) %>%
  kibble_one(format = 'html') %>%
  kableExtra::kable_styling(bootstrap_options = 'striped') %>% 
  print()







