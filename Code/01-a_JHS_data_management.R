
rm(list=ls())

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

# echa.sas7bdat
# lvedd, pwt, swt 
# for m-mode, ECHA43, ECHA45, ECHA41  
# for 2D,     ECHA52, ECHA54, ECHA50

echo<-haven::read_sas("Datasets/JHS/echa.sas7bdat")%>%
  set_names(tolower(names(.)))%>%
  mutate(lvm_2d=0.8*(1.04*(((echa50/10)+(echa52/10)+(echa54/10))^3-(echa52/10)^3))+0.6,
         lvm_mm=0.8*(1.04*(((echa41/10)+(echa43/10)+(echa45/10))^3-(echa43/10)^3))+0.6)%>%
  dplyr::select(subjid,lvm_2d)%>%
  na.omit()

library(tidyverse)
library(magrittr)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) > 0
}

anly = read.csv('Datasets/JHS/JHS_baseline_visit.csv',
                stringsAsFactors=FALSE)%>%
  mutate(bsa=sqrt((height*weight)/3600))%>%
  dplyr::filter(subjid%in%echo$subjid)%>%
  merge(echo,by='subjid')%>%
  mutate(lvm=lvm_2d/bsa,
         lvh=ifelse((sex=='Male' & lvm > 115) | 
                    (sex=='Female' & lvm > 95), 1, 0),
         race=1, study=1) %>%
  dplyr::select(subjid,lvm,lvh,sex,cln_sbp=sbp,cln_dbp=dbp,
                age,alc,currentSmoker,diab=Diabetes,
                edu=edu3cat,race,study)

for(i in names(anly)){
  if(is.character(anly[[i]])){
    anly[[i]][anly[[i]]%in%c('',' ')]=NA
    anly[[i]]%<>%factor()
  } 
}

levels(anly$edu)<-list(
  '1'="High school graduate/GED",
  '1'="Less than high school",
  '2'="Attended vocational school, trade school, or college"
)

levels(anly$currentSmoker)<-list('1'='No','2'='Yes')
levels(anly$diab)<-list('1'='No','2'='Yes')
levels(anly$sex)<-list('1'='Female','2'='Male')

anly%<>%mutate(alc=as.numeric(alc),
               currentSmoker=as.numeric(currentSmoker),
               diab=as.numeric(diab),
               edu=as.numeric(edu),
               sex=as.numeric(sex))


labs=nobs=list()
labs[[1]]="JHS participants"
nobs[[1]]=nrow(anly)
exc=y30_exc=list(labs=labs,nobs=nobs)

exc$labs[[2]]="Participants who underwent 24-hour ABPM"
exc$nobs[[2]]=1148

exc$labs[[3]]<-
  "Participants with \u2265 10 awake and \u2265 5 asleep BP measurements."
exc$nobs[[3]]=1046

exc$labs[[4]]="Participants who provided asleep and awake times"
exc$nobs[[4]]=943

abpm = read.csv("Datasets/JHS/JHS_abpm_943_long.csv") %>% 
  dplyr::select(-tsw,-hour,-minute,-nt,-ends_with('_cal'),
                -time_seconds_since_midnight)%>%
  mutate(nht=ifelse(slp_sbp>=120 | slp_dbp>=70,1,0),
         slp_dur=NA) %>%
  mutate_cond(time_awake<=time_asleep&!is.na(time_awake)&!is.na(time_asleep), 
              slp_dur=(24-time_asleep)+time_awake) %>%
  mutate_cond(time_awake>time_asleep&!is.na(time_awake)&!is.na(time_asleep), 
              slp_dur=time_awake-time_asleep)%>%
  left_join(anly,by='subjid')%>%
  plyr::ddply(.variables = 'subjid',.fun=function(blk){
    if(any(blk$awake[blk$time%>%is.between(1,5)]==1)){
      #cat(blk$subjid[1]); 
      return(NULL)
    } else {
      return(blk)
    }
  })%>% 
  dplyr::rename(tsm=time)%>%
  dplyr::select(-ppr,-hr,-slp_state)

exc$labs[[5]]="Participants who were asleep during sampling times (1am-5am)"
exc$nobs[[5]]=length(unique(abpm$subjid))

abpm%<>%mutate_if(is.factor,as.numeric)

library(lme4)
library(splines)
set.seed(329)

lmm<-lme4::lmer(sbp~bs(tsm,4)+(bs(tsm,4)||subjid),data=abpm)
abpm$sbp_imp<-predict(lmm)+rnorm(nrow(abpm),mean=0,sd=sigma(lmm))

lmm<-lme4::lmer(dbp~bs(tsm,4)+(bs(tsm,4)||subjid),data=abpm)
abpm$dbp_imp<-predict(lmm)+rnorm(nrow(abpm),mean=0,sd=sigma(lmm))

abpm%<>%group_by(subjid)%>%
  mutate(slp_sbp_imp=mean(sbp_imp[awake==0]),
         slp_dbp_imp=mean(dbp_imp[awake==0]))

abpm_smry<-abpm%>%ungroup()%>%
  dplyr::select(subjid,slp_sbp,slp_sbp_imp,slp_dbp,slp_dbp_imp,awk_sbp,
                awk_dbp,slp_dur,nht)%>%unique()%>%
  mutate(subjid=as.character(subjid))

anly%<>%dplyr::filter(subjid%in%abpm$subjid)%>%
  merge(abpm_smry,by='subjid')%>%
  mutate(subjid=as.character(subjid))

saveRDS(abpm,'Datasets/JHS_abpm.RDS')
saveRDS(anly,'Datasets/JHS_anly.RDS')
saveRDS(exc ,'Datasets/JHS_excl.RDS')






