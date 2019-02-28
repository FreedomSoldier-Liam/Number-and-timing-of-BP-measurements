

library(tidyverse)
library(magrittr)
library(haven)
library(splines)

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
  "jhs_abpm_long_27FEB2019.csv") %>%
  read_csv() %>%
  dplyr::filter(sr_journal_provided == 1) %>% # 1015
  dplyr::select(
    -ppr,
    -hr, 
    -sr_journal_provided,
    -tsw, 
    -dt, 
    -nt
  ) %>%
  group_by(subjid) %>% 
  mutate(
    asleep_1_to_5 = all(awake[time>=1 & time<=5]==0),
    nreadings_awk = sum(awake==1),
    nreadings_slp = sum(awake==0),
    slp_sbp = mean(sbp[awake==0]),
    slp_dbp = mean(dbp[awake==0])
  ) 

exc$labs[[3]] = "Participants who reported times of falling asleep and waking"
exc$nobs[[3]] = length(unique(jhs_abpm_long$subjid))

jhs_abpm_long %<>% 
  dplyr::filter(
    nreadings_awk >=10,
    nreadings_slp >= 5
  )

exc$labs[[4]] <-
  "Participants with \u2265 10 awake and \u2265 5 asleep BP measurements."
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


# Pseudo asleep BP means --------------------------------------------------


lmm <- lme4::lmer(
  sbp~bs(time,4)+(bs(time,4)||subjid),
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
  dbp~bs(time,4)+(bs(time,4)||subjid),
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
    alc,
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
  "jhs_abpm_means_27FEB2019.csv"
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
)

saveRDS(jhs_abpm_long,'Datasets/JHS_abpm_long.RDS')
saveRDS(analysis,'Datasets/JHS_analysis.RDS')
saveRDS(exc ,'Datasets/JHS_excl.RDS')







