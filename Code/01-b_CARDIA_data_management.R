rm(list=ls())

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
library(tidyverse)
library(haven)
library(lubridate)
library(magrittr)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) > 0
}

echo<-read_sas("Datasets/CARDIA/cardia_echo.sas7bdat")%>%
  dplyr::select(id=ID,lvmass2d=ILVMASS2D)%>%
  dplyr::mutate(id=substr(id,0,5),id=as.numeric(id))

demo<- haven::read_sas('Datasets/CARDIA/cardia_abpm_781_updated.sas7bdat')%>%
  set_names(tolower(names(.)))%>%
  mutate(edu=cut(educ,breaks = c(0,12,Inf),labels=c(1,2)),
         edu=as.numeric(edu),
         diabetes=diabetes+1,
         currentSmoker=ifelse(smoke==0,2,1),
         alc=ifelse(alcohol==2,1,2))%>%
  dplyr::select(id=short_id, age, sex=sex_cat, race=race_cat, 
                height_cm, weight_lbs,edu,diab=diabetes,alc,
                currentSmoker,edu,cln_sbp=c_sbp,cln_dbp=c_dbp)%>%
  data.frame()

# location to read in data
schwartz_abpm_loc=file.path('O:','REGARDS','CARDIA','Data',
                            'Population','ABPM_Schwartz')

uab_slp=read_sas(file.path(
  schwartz_abpm_loc,'final_summary_dataset_uab_v1.sas7bdat'))
nw_slp<-read_sas(file.path(
  schwartz_abpm_loc,'final_summary_dataset_nw_v1.sas7bdat'))

slp=rbind(uab_slp,nw_slp)%>%
  mutate(slp_dur=as.numeric(sleep_duration)/60^2)%>%
  dplyr::select(id,slp_dur)

demo%<>%
  left_join(slp,by='id')%>%
  left_join(echo,by='id')%>%
  mutate(weight_kg=weight_lbs/2.20462,
         bsa=sqrt((height_cm*weight_kg)/3600),
         lvm=lvmass2d/bsa,
         lvh=ifelse((sex==2 & lvm > 115) | 
                      (sex==1 & lvm > 95), 1, 0),
         study=2) %>%
  dplyr::select(-lvmass2d,-height_cm,-weight_lbs,-weight_kg,-bsa)

demo$sex[demo$sex==0]=2
demo$race[demo$race==0]=2

abpm<-read_sas("Datasets/CARDIA/spacelabs_abp2_edited.sas7bdat")%>%
  set_names(tolower(names(.)))%>%
  data.frame()

abpm$period[abpm$period=='WU']='SL'

suppressWarnings(
  abpm%<>%dplyr::filter(period%in%c("AW","SL"),
                        valid_reading==1,
                        manual_reading==0)%>%
    dplyr::select(id,time=tim,period,sbp=sysbp,dbp=diabp)%>%
    mutate(time=as.numeric(hms(time))/(60^2))%>%
    group_by(id)%>%
    mutate(nawake=sum(period=="AW"),
           nasleep=sum(period=='SL'),
           awk_sbp=mean(sbp[period=='AW']),
           slp_sbp=mean(sbp[period=='SL']),
           awk_dbp=mean(dbp[period=='AW']),
           slp_dbp=mean(dbp[period=='SL']),
           nht=ifelse(slp_sbp>=120 | slp_dbp>=70,1,0),
           period=factor(period))%>%
    dplyr::filter(nawake>10,nasleep>5)%>%
    group_by(id,period)%>%
    mutate(tss=time-min(time))%>%
    mutate_cond(period=="AW",tss=0)%>%
    ungroup()%>%
    mutate(id=as.character(id),
           awake=ifelse(period=='AW',1,0))%>%  
    plyr::ddply(.variables = 'id',.fun=function(blk){
      if(any(blk$awake[blk$time%>%is.between(1,5)]==1)){
        return(NULL)
      } else {
        return(blk)
      }
    })
)

abpm$tss[abpm$awake==1]=NA

#View(abpm[abpm$id==abpm$id[which.max(abpm$tss)],])

abpm$time[abpm$time>=24]%<>%subtract(24)
abpm$time[abpm$time>=24]%<>%subtract(24)

ids_to_filter=abpm$id[(abpm$period=='AW'&is.between(abpm$time,1,5))]
ids_to_filter_counts=table(ids_to_filter)
ids_to_filter=names(ids_to_filter_counts)[ids_to_filter_counts>1]


abpm%<>%dplyr::filter(!(id%in%ids_to_filter))%>%
  mutate(id=as.numeric(id))%>%
  merge(demo,by='id')%>%
  dplyr::rename(subjid=id,tsm=time)

set.seed(329)

lmm<-lme4::lmer(sbp~bs(tsm,4)+(bs(tsm,4)||subjid),data=abpm)
abpm$sbp_imp<-predict(lmm)+rnorm(nrow(abpm),mean=0,sd=sigma(lmm))

lmm<-lme4::lmer(dbp~bs(tsm,4)+(bs(tsm,4)||subjid),data=abpm)
abpm$dbp_imp<-predict(lmm)+rnorm(nrow(abpm),mean=0,sd=sigma(lmm))

abpm%<>%group_by(subjid)%>%
  mutate(slp_sbp_imp=mean(sbp_imp[awake==0]),
         slp_dbp_imp=mean(dbp_imp[awake==0]))


labs=nobs=list()
labs[[1]]="CARDIA participants"
nobs[[1]]=5115
exc=y30_exc=list(labs=labs,nobs=nobs)

exc$labs[[2]]="Participants who underwent 24-hour ABPM"
exc$nobs[[2]]=825

exc$labs[[3]]<-
  "Participants with \u2265 10 awake and \u2265 5 asleep BP measurements."
exc$nobs[[3]]=781

exc$labs[[4]]="Participants who provided asleep and awake times"
exc$nobs[[4]]=781

exc$labs[[5]]="Participants who were asleep during sampling times (1am-5am)"
exc$nobs[[5]]=length(unique(abpm$subjid))

abpm%<>%dplyr::select(subjid,tsm,tss,awake,sbp,dbp,slp_sbp,slp_dbp,
                      awk_sbp,awk_dbp,nht,lvm,lvh,sex,age,alc,
                      slp_sbp_imp,slp_dbp_imp,sbp_imp,dbp_imp,
                      currentSmoker,diab,edu,race,study)

mvars=c("subjid","slp_sbp","slp_dbp","awk_sbp","awk_dbp","nht",
        "slp_sbp_imp","slp_dbp_imp")

demo%<>%dplyr::rename(subjid=id)%>%
  mutate(subjid=as.numeric(subjid))%>%
  right_join(unique(abpm[,mvars]),by='subjid')%>%
  dplyr::select(subjid,lvm,lvh,sex,age,alc,currentSmoker,slp_dur,
                diab,edu,race,study,slp_sbp,slp_dbp,
                slp_sbp_imp,slp_dbp_imp,
                awk_sbp,awk_dbp,nht,cln_sbp,cln_dbp)

saveRDS(abpm,'Datasets/CARDIA_abpm.RDS')
saveRDS(demo,'Datasets/CARDIA_anly.RDS')
saveRDS(exc,'Datasets/CARDIA_exc.RDS')

# write.csv(abpm,"../Change in ABPM over 25 years/Clean Data/CARDIA_abpm_y30.csv")
## Quality control
# jo_means=read_sas("Datasets/CARDIA/y30_abpm_means.sas7bdat")%>%
#   dplyr::filter(GE_5_Valid_sleep==1,GE_10_Valid_awake==1)%>%
#   dplyr::select(id,mn_sysbp_AW,mn_sysbp_SL,mn_diabp_AW,mn_diabp_SL,
#                 time_asleep=Sleep_Time,time_awake=Wakeup_Time)%>%
#   mutate(time_asleep=as.numeric(hms(time_asleep))/(60^2),
#          time_awake=as.numeric(hms(time_awake))/(60^2))
# 
# my_means=unique(abpm[,c("id","awk_sbp","slp_sbp")])
# 
# (my_means$slp_sbp-jo_means$mn_sysbp_SL)%>%sort(decreasing = T)%>%
#   magrittr::extract(1)
# 
# id_for_inspect=jo_means$id[which.max(my_means$slp_sbp-jo_means$mn_sysbp_SL)]
# 
# abpm[abpm$id==id_for_inspect,]%>%data.frame()
