
library(tidyverse)
library(haven)
library(lubridate)
library(magrittr)


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



labs=nobs = list()
labs[[1]] = "CARDIA participants"
nobs[[1]] = 5115
exc       = list(labs=labs,nobs=nobs)

cardia_abpm_means <- file.path(
  data_file_path,
  "CARDIA ABPM",
  "Year 30",
  "CARDIA_ABPM_Y30_means.csv"
) %>% 
  read_csv() %>% 
  dplyr::filter(
    nreadings_asleep >= 5
  )

cardia_abpm_long <- file.path(
  data_file_path,
  "CARDIA ABPM",
  "Year 30",
  "CARDIA_ABPM_Y30_long.csv"
) %>% 
  read_csv() %>% 
  set_names(
    tolower(names(.))
  ) %>% 
  dplyr::select(
    -c(
      center,
      order,
      hr,
      daytime,
      nighttime
    )
  ) %>% 
  dplyr::rename(tsm = time) %>% 
  dplyr::arrange(sid, tsm)

exc$labs[[2]]="Participants who underwent 24-hour ABPM"
exc$nobs[[2]]=length(unique(cardia_abpm_long$sid))

cardia_abpm_long %<>% 
  dplyr::filter(sid %in% cardia_abpm_means$sid) 

exc$labs[[3]]<-
  "Participants with \u2265 5 asleep BP measurements."
exc$nobs[[3]]=length(unique(cardia_abpm_long$sid))

cardia_abpm_long %<>% 
  group_by(sid) %>% 
  mutate(
    asleep_1_to_5 = all(period[tsm>=1 & tsm<=5]%in%c("SL","WU"))
  ) %>%
  dplyr::filter(
    asleep_1_to_5 == TRUE,
    period %in% c("AW","SL")
  ) %>%
  group_by(sid, period) %>% 
  mutate(
    tss = tsm - tsm[1]
  ) %>% 
  ungroup() %>% 
  mutate(
    tss = case_when(
      period == 'AW' ~ 0,
      tss < 0 ~ tss + 24,
      TRUE ~ tss
    )
  )

exc$labs[[4]]="Participants who were asleep during sampling times (1am-5am)"
exc$nobs[[4]]=length(unique(cardia_abpm_long$sid))

valid_samples <- cardia_abpm_long %>% 
  group_by(sid) %>% 
  nest() %>% 
  mutate(
    abp_tsm_valid = map_lgl(
      data, 
      .f=abp_sampler, 
      times=seq(from = 1, to = 5, by = 1/2), 
      time_variable = 'tsm', 
      quality_check=TRUE
    ),
    abp_tss_valid = map_lgl(
      data, 
      .f = abp_sampler, 
      times = seq(from = 1, to = 5, by = 1/2),
      time_variable = 'tss',
      quality_check = TRUE
    )
  ) %>% 
  dplyr::filter(
    abp_tsm_valid==TRUE,
    abp_tss_valid==TRUE
  ) %>% 
  pluck('sid')

cardia_abpm_long %<>% 
  dplyr::filter(
    sid %in% valid_samples
  )

exc$labs[[5]]="Participants who had a BP reading within 30 minutes of all sampling times"
exc$nobs[[5]]=length(unique(cardia_abpm_long$sid))

set.seed(329)

lmm<-lme4::lmer(
  sbp~bs(tsm)+(bs(tsm)|sid),
  data=cardia_abpm_long
)

cardia_abpm_long$sbp_imp <- predict(lmm) %>% 
  add(
    rnorm(
      length(.),
      mean=0,
      sd=sigma(lmm)
    )
  )

lmm <- lme4::lmer(
  dbp~bs(tsm)+(bs(tsm)|sid),
  data=cardia_abpm_long
)

cardia_abpm_long$dbp_imp <- predict(lmm) %>% 
  add(
    rnorm(
      length(.),
      mean=0,
      sd=sigma(lmm)
    )
  )

cardia_abpm_long %<>% 
  dplyr::select(
    subjid = sid, 
    sbp,
    sbp_imp,
    dbp, 
    dbp_imp,
    tss,
    tsm
  )

echo <- file.path(
  data_file_path, 
  'CARDIA_analysis',
  'CARDIA_ECHO_Y30.csv'
) %>%
  read_csv() 

cardia_y30 <- file.path(
  data_file_path,
  "CARDIA_analysis",
  "cardia_analysis_all_visits.csv"
) %>%
  read_csv() %>% 
  dplyr::filter(
    exam == 30,
    sid %in% cardia_abpm_long$subjid
  ) %>% 
  left_join(echo, by = 'sid') %>% 
  left_join(cardia_abpm_means, by = 'sid') %>% 
  dplyr::select(
    subjid = sid,
    age,
    sex,
    race,
    edu = educ,
    diabetes=diab,
    currentsmoker = smoke,
    bpmeds = bpmeds,
    cln_sbp = sbp,
    cln_dbp = dbp,
    lvm = lvm_bsa,
    lvh = lvh_bsa,
    acr,
    slp_sbp,
    slp_dbp,
    awk_sbp,
    awk_dbp,
    sleep_duration
  ) %>% 
  mutate(
    nht=ifelse(slp_sbp>=120 | slp_dbp>=70,1,0),
    currentsmoker = case_when(
      currentsmoker %in% c('Never','Former') ~ "No",
      currentsmoker =='Current' ~ "Yes"
    ),
    subjid = as.character(subjid)
  ) %>% 
  dplyr::filter(
    sleep_duration >= 5
  )

exc$labs[[6]]="Participants who slept for \u2265 5 hours"
exc$nobs[[6]]=nrow(cardia_y30)

cardia_abpm_long %<>% 
  ungroup() %>% 
  mutate(
    subjid=as.character(subjid)
  ) %>% 
  dplyr::filter(
    subjid %in% cardia_y30$subjid
  ) %>% 
  left_join(
    cardia_y30[,c("subjid","race","sex")],
    by = 'subjid'
  )

saveRDS(cardia_abpm_long,'Datasets/CARDIA_abpm_long.RDS')
saveRDS(cardia_y30,'Datasets/CARDIA_analysis.RDS')
saveRDS(exc ,'Datasets/CARDIA_excl.RDS')
