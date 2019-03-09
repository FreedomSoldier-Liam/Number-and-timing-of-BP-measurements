
library(kableExtra)
library(tidyverse)
library(magrittr)

files <- list.files(path='Code/KableWon',all.files=TRUE, 
                    full.names=TRUE, pattern='.R')

for(f in files) source(f)

real_data=bind_rows(
  jhs    = readRDS('Datasets/JHS_analysis.RDS'),
  cardia = readRDS('Datasets/CARDIA_analysis.RDS'),
  .id = "study"
) %>% 
  mutate(
    edu = case_when(
      edu %in% c(
        'Attended vocational school, trade school, or college',
        'College',
        'Graduate School'
      ) ~ "College or professional degree",
      edu %in% c(
        "High school graduate/GED",
        "High School"
      ) ~ "High School graduate",
      edu %in% c(
        'Less than high school',
        'Jr. High'
      )  ~ "Did not graduate from high school"
    ),
    age = case_when(
      study=='cardia' ~ age + 30,
      study=='jhs' ~ age
    )
  ) %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(
    edu = factor(
      edu, 
      levels = c(
        "College or professional degree", 
        "High School graduate",
        "Did not graduate from high school"
      )
    ),
    nht = factor(
      nht, 
      levels = 0:1,
      labels = c(
        "Nocturnally normotensive",
        "Nocturnally hypertensive"
      )
    ),
    lvh = factor(
      lvh, 
      levels = 0:1,
      labels = c(
        'No ventricular hypertrophy',
        'Left ventricular hypertrophy'
      )
    ),
    study = fct_relevel(study, "jhs")
  ) 
  

meta_data=read.csv('Datasets/meta_data.csv',stringsAsFactors=F)

meta = real_data %>% 
  mutate_if(is.factor, as.numeric) %>% 
  as.data.frame() %>% 
  create_meta_object(
    meta_data = meta_data,
    exposure  = 'study'
  )

tab1=make_table_one(meta=meta,as.kable=TRUE,
                    include.pval = FALSE)

print(tab1)