

library(tidyverse)
library(magrittr)
library(haven)
library(splines)

files <- list.files(
  path='Code/abp_sampling_functions/',
  all.files=TRUE, 
  full.names=TRUE, pattern='.R'
)

for(f in files) source(f)

data_file_path <- file.path(
  "..",
  "..",
  "..",
  "..",
  "REGARDS",
  "JHS",
  "Derived datasets",
  "04-06-2019",
  "data",
  "output"
)

labs=nobs=list()
labs[[1]]="JHS participants"
nobs[[1]]=5306
exc=list(labs=labs,nobs=nobs)

exc$labs[[2]]="Participants who underwent 24-hour ABPM"
exc$nobs[[2]]=1148



jhs_abpm_long <- file.path(
  data_file_path,
  "abpm_mrp.csv"
) %>%
  read_csv() %>%
  dplyr::filter(sleep_diary_valid == 1) %>% 
  group_by(subjid) %>% 
  mutate(
    tss = time_msr - time_slp,
    tss = case_when(
      status!='Asleep' ~ 0,
      tss < 0 ~ tss + 24,
      TRUE ~ tss
    )
  ) %>% 
  dplyr::rename(tsm = time_msr) %>% 
  group_by(subjid) %>% 
  mutate(
    asleep_1_to_5 = all(status[tsm>=1 & tsm<=5]=="Asleep"),
    nreadings_awk = sum(status=='Awake'),
    nreadings_slp = sum(status=='Asleep')
    #max_diff_tsm = max(diff(tss))
  )

exc$labs[[3]] = "Participants with valid sleep diaries"
exc$nobs[[3]] = length(unique(jhs_abpm_long$subjid))

jhs_abpm_long %<>% 
  dplyr::filter(nreadings_slp >= 5)

exc$labs[[4]] <-"Participants with \u2265 5 asleep BP measurements."
exc$nobs[[4]] = length(unique(jhs_abpm_long$subjid))

jhs_abpm_long %<>% 
  dplyr::filter(asleep_1_to_5 == 1) %>% 
  dplyr::select(
    -starts_with("nread"), 
    -asleep_1_to_5
  )

exc$labs[[5]]="Participants who were asleep during sampling times (1am-5am)"
exc$nobs[[5]]=length(unique(jhs_abpm_long$subjid))

valid_samples <- jhs_abpm_long %>% 
  group_by(subjid) %>% 
  nest() %>% 
  mutate(
    abp_tsm_valid = map_lgl(
      data, 
      .f=abp_sampler, 
      times=seq(from = 1, to = 5, by = 1/3), 
      time_variable = 'tsm', 
      quality_check=TRUE
    ),
    abp_tss_valid = map_lgl(
      data, 
      .f = abp_sampler, 
      times = seq(from = 1, to = 5, by = 1/3),
      time_variable = 'tss',
      quality_check = TRUE
    )
  ) %>% 
  dplyr::filter(
    abp_tsm_valid==TRUE,
    abp_tss_valid==TRUE
  ) %>% 
  pluck('subjid')

jhs_abpm_long %<>% 
  dplyr::filter(subjid %in% valid_samples)

exc$labs[[6]]="Participants who had a BP reading within 30 minutes of all sampling times"
exc$nobs[[6]]=length(unique(jhs_abpm_long$subjid))


# Pseudo asleep BP means --------------------------------------------------

nimpute = 10
set.seed(329)

lmm <- lme4::lmer(
  sbp~bs(tsm)+(bs(tsm)|subjid),
  data=jhs_abpm_long
)

for(i in 1:nimpute){
  jhs_abpm_long[[paste0("sbp_imp_",i)]] <- 
    predict(lmm) %>% 
    add(
      rnorm(
        length(.),
        mean=0,
        sd=sigma(lmm)
      )
    )
}

lmm <- lme4::lmer(
  dbp~bs(tsm)+(bs(tsm)|subjid),
  data=jhs_abpm_long
)

for(i in 1:nimpute){
  jhs_abpm_long[[paste0("dbp_imp_",i)]] <- 
    predict(lmm) %>% 
    add(
      rnorm(
        length(.),
        mean=0,
        sd=sigma(lmm)
      )
    )
}

jhs_abpm_long %<>% 
  dplyr::select(
    subjid, 
    tsm, 
    tss,
    sbp,
    dbp, 
    starts_with('sbp_imp'),
    starts_with('dbp_imp')
  )

is_01 <- function(x){
  all(na.omit(x) %in% c(0,1))
}

jhs_vis1 <- file.path(
  "..","..","Datasets",
  "JHS_analysis",
  "Processed data",
  "jhs_visit1.csv"
) %>%
  read_csv() %>% 
  dplyr::filter(
    subjid %in% jhs_abpm_long$subjid
  ) %>% 
  mutate(race='Black') %>% 
  dplyr::select(
    subjid,
    age,
    sex,
    race,
    edu=edu3cat,
    diabetes,
    #alc,
    currentsmoker,
    bpmeds,
    cln_sbp=sbp,
    cln_dbp=dbp,
    lvm,
    lvh,
    acr=albumin_creatinine_ratio
  ) %>% 
  dplyr::mutate(
    albuminuria = ifelse(acr > 30, 1, 0),
    edu = factor(edu,
      levels = c(0,1,2),
      labels = c(
        "Less than High School",
        "High School graduate/GED",
        "College graduate"
      ) 
    )
  )

jhs_abpm_means <- file.path(
  data_file_path,
  "abpm_srp.csv"
) %>%
  read_csv() %>%
  dplyr::filter(subjid %in% jhs_abpm_long$subjid) %>% 
  dplyr::select(
    subjid, 
    slp_sbp,
    slp_dbp,
    awk_sbp,
    awk_dbp,
    slp_duration
  ) %>% 
  group_by(subjid) %>%
  mutate(nht=ifelse(slp_sbp>=120 | slp_dbp>=70,1,0))

analysis <- dplyr::left_join(
  jhs_vis1,
  jhs_abpm_means,
  by='subjid'
) %>% 
  dplyr::filter(
    slp_duration >= 5 & slp_duration < 12
  ) %>% 
  mutate_if(is_01, factor, labels = c("No","Yes")) %>% 
  mutate_if(is.factor, as.character)

exc$labs[[7]]="Participants who slept for \u2265 5 and < 12 hours"
exc$nobs[[7]]=nrow(analysis)

jhs_abpm_long %<>% 
  ungroup() %>% 
  mutate(subjid=as.character(subjid)) %>% 
  dplyr::filter(subjid %in% analysis$subjid) %>% 
  left_join(
    analysis[,c("subjid","race","sex")],
    by = 'subjid'
  )

saveRDS(jhs_abpm_long,'Datasets/JHS_abpm_long.RDS')
saveRDS(analysis,'Datasets/JHS_cleaned.RDS')
saveRDS(exc ,'Datasets/JHS_excl.RDS')







