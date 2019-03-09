

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
  "Datasets"
)

labs=nobs=list()
labs[[1]]="JHS participants"
nobs[[1]]=5306
exc=list(labs=labs,nobs=nobs)

exc$labs[[2]]="Participants who underwent 24-hour ABPM"
exc$nobs[[2]]=1148


jhs_abpm_long <- file.path(
  data_file_path,
  "JHS ABPM",
  "long",
  "jhs_abpm_long_27FEB2019.csv"
) %>%
  read_csv() %>%
  dplyr::filter(
    sr_journal_provided == 1
  ) %>% # 1015
  dplyr::select(
    -ppr,
    -hr, 
    -starts_with("sr_"),
    -tsw, 
    -dt, 
    -nt
  ) %>%
  dplyr::rename(
    tsm = time
  ) %>% 
  mutate(
    tss = case_when(
      awake==1 ~ 0,
      awake==0 ~ tss
    )
  ) %>% 
  group_by(subjid) %>% 
  mutate(
    asleep_1_to_5 = all(awake[tsm>=1 & tsm<=5]==0),
    nreadings_awk = sum(awake==1),
    nreadings_slp = sum(awake==0),
    slp_sbp = mean(sbp[awake==0]),
    slp_dbp = mean(dbp[awake==0])
    #max_diff_tsm = max(diff(tss))
  )

exc$labs[[3]] = "Participants who reported times of falling asleep and waking"
exc$nobs[[3]] = length(unique(jhs_abpm_long$subjid))

jhs_abpm_long %<>% 
  dplyr::filter(
    nreadings_slp >= 5
  )

exc$labs[[4]] <-
  "Participants with \u2265 5 asleep BP measurements."
exc$nobs[[4]] = length(unique(jhs_abpm_long$subjid))

jhs_abpm_long %<>% 
  dplyr::filter(
    asleep_1_to_5 == 1
  ) %>% 
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
  dplyr::filter(
    subjid %in% valid_samples
  )

exc$labs[[6]]="Participants who had a BP reading within 30 minutes of all sampling times"
exc$nobs[[6]]=length(unique(jhs_abpm_long$subjid))


# Pseudo asleep BP means --------------------------------------------------

lmm <- lme4::lmer(
  sbp~bs(tsm)+(bs(tsm)|subjid),
  data=jhs_abpm_long
)

jhs_abpm_long$sbp_imp <- predict(lmm) %>% 
  add(
    rnorm(
      length(.),
      mean=0,
      sd=sigma(lmm)
    )
  )

lmm <- lme4::lmer(
  dbp~bs(tsm)+(bs(tsm)|subjid),
  data=jhs_abpm_long
)

jhs_abpm_long$dbp_imp <- predict(lmm) %>% 
  add(
    rnorm(
      length(.),
      mean=0,
      sd=sigma(lmm)
    )
  )

jhs_abpm_long %<>% 
  dplyr::select(
    subjid, 
    sbp,
    sbp_imp,
    dbp, 
    dbp_imp,
    tsm, 
    tss
  )

jhs_vis1 <- file.path(
  data_file_path,
  "JHS_analysis",
  "Processed data",
  "jhs_visit1.csv") %>%
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
  ) 

jhs_abpm_means <- file.path(
  data_file_path,
  "JHS ABPM",
  "means",
  "jhs_abpm_means_05MAR2019.csv"
) %>%
  read_csv() %>%
  dplyr::filter(subjid %in% jhs_abpm_long$subjid) %>% 
  dplyr::select(
    -c(
      starts_with("dt"),
      starts_with("nt"),
      starts_with("nread")
    )
  ) %>% 
  group_by(subjid) %>%
  mutate(nht=ifelse(slp_sbp>=120 | slp_dbp>=70,1,0))

analysis <- dplyr::left_join(
  jhs_vis1,
  jhs_abpm_means,
  by='subjid'
) %>% 
  dplyr::filter(
    sleep_duration >= 5
  )

exc$labs[[7]]="Participants who slept for \u2265 5 hours"
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
saveRDS(analysis,'Datasets/JHS_analysis.RDS')
saveRDS(exc ,'Datasets/JHS_excl.RDS')







