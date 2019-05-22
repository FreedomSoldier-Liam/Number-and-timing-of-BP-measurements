library(kableExtra)
library(tidyverse)
library(magrittr)
library(KableOne)
library(labelled)


data <- bind_rows(
  jhs    = readRDS('Datasets/JHS_cleaned.RDS'),
  cardia = readRDS('Datasets/CARDIA_cleaned.RDS'),
  .id = "study"
) %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(
    age = case_when(
      study=='cardia' ~ age + 30,
      study=='jhs' ~ age
    ),
    study = fct_recode(study, 
      'JHS' = 'jhs',
      'CARDIA' = 'cardia'
    ),
    study = fct_relevel(study, "JHS")
  ) %>%
  select_labelled(
    subjid = 'Participant identification',
    age	= 'Age',
    sex	= 'Sex',
    race = 'Race',
    edu = 'Education',
    currentsmoker =	'Current smoker',
    diabetes = 'Diabetes',
    albuminuria = 'Albuminuria',
    lvm	= 'Left ventricular mass indexed to BSA',
    lvh = 'Left ventricular hypertrophy',
    slp_duration = 'Sleep duration',
    nht = 'Nocturnal hypertension',
    bpmeds = 'Antihypertensive medication use',
    slp_sbp =	'Nighttime systolic',
    slp_dbp	= 'Nighttime diastolic',
    awk_sbp =	'Daytime systolic',
    awk_dbp =	'Daytime diastolic',
    cln_sbp =	'Clinic systolic',
    cln_dbp =	'Clinic diastolic',
    study	= 'Study'
  ) %>% 
  set_variable_groups(
    `Blood pressure, mm Hg` = c(
      'slp_sbp','slp_dbp',
      'awk_sbp','awk_dbp',
      'cln_sbp','cln_dbp'
    )
  ) %>% 
  set_variable_abbrs(
    lvm = "BSA = body surface area",
    study = c(
      "JHS = Jackson Heart Study",
      "CARDIA = Coronary Artery Risk Development in Young Adults"
    )
  ) %>% 
  set_variable_units(
    age = 'years',
    slp_duration = 'hours',
    lvm = 'g/m<sup>2</sup>'
  ) %>% 
  set_variable_notes(
    albuminuria = 'Albuminuria was defined as urinary albumin-to-creatinine ratio  > 30 mg/g'
  )

write_rds(
  data,
  "Datasets/pooled_analysis.RDS"
)

write_rds(
  filter(data, study == 'JHS'),
  "Datasets/JHS_analysis.RDS"
)

write_rds(
  filter(data, study == 'CARDIA'),
  "Datasets/CARDIA_analysis.RDS"
)