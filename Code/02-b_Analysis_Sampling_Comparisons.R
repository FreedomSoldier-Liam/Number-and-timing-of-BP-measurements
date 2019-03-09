mae<-function(x,y) mean(abs(x-y),na.rm=TRUE)

library(tidyverse)
library(magrittr)
library(fmsb)
library(BlandAltmanLeh)

files <- list.files(path='Code/abp_sampling_functions',all.files=TRUE, 
  full.names=TRUE, pattern='.R')

for(f in files) source(f)


## Distributed protocol times
lwr_time=1
upr_time=5
ntimes=list(2,3,4)

distr = map(
  ntimes, ~combn(lwr_time:upr_time,.)
) %>%
  map(
    .f=function(mat){
      map(1:ncol(mat),~mat[,.]
      )
    }
  ) %>%
  set_names(
    paste("msr",ntimes,'times',sep='_')
  )

cnctr=list(
  JHS=list(
    rep(1,4)+c(0,.33,.66, 1),
    rep(2,4)+c(0,.33,.66, 1),
    rep(3,4)+c(0,.33,.66, 1),
    rep(4,4)+c(0,.33,.66, 1),
    rep(1,3)+c(0,.33,.66),
    rep(2,3)+c(0,.33,.66),
    rep(3,3)+c(0,.33,.66),
    rep(4,3)+c(0,.33,.66),
    rep(1,2)+c(0,.33),
    rep(2,2)+c(0,.33),
    rep(3,2)+c(0,.33),
    rep(4,2)+c(0,.33)
  ),
  CARDIA=list(
    rep(1,4)+c(0,1/2,1,3/2),
    rep(2,4)+c(0,1/2,1,3/2),
    rep(3,4)+c(0,1/2,1,3/2),
    rep(4,4)+c(0,1/2,1,3/2),
    rep(1,3)+c(0,1/2,1),
    rep(2,3)+c(0,1/2,1),
    rep(3,3)+c(0,1/2,1),
    rep(4,3)+c(0,1/2,1),
    rep(1,2)+c(0,1/2),
    rep(2,2)+c(0,1/2),
    rep(3,2)+c(0,1/2),
    rep(4,2)+c(0,1/2)
  )
)

inputs=expand.grid(
  data_lab=c("JHS","CARDIA"),
  impute=c(FALSE,TRUE),
  strata=c("overall",'sex','race'),
  tvarbl=c('tss','tsm'),
  stringsAsFactors = FALSE
)

take_out=which(inputs$strata=='race' & inputs$data_lab=='JHS')
inputs=inputs[-take_out,]

means_data <- 
  list(
    JHS=readRDS("Datasets/JHS_analysis.RDS"),
    CARDIA=readRDS("Datasets/CARDIA_analysis.RDS")
  ) %>% 
  map(
    ~{
      dplyr::select(.x, subjid, starts_with("slp_"), nht) %>% 
        dplyr::mutate(subjid = as.character(subjid))
    }
  )

for(i in 1:nrow(inputs)){
  
  input=inputs[i,]
  
  abpm_long <- readRDS(
    paste0("Datasets/",input$data_lab,"_abpm_long.RDS")
  )
  
  if(input$impute){
    abpm_long%<>%dplyr::mutate(sbp=sbp_imp,dbp=dbp_imp)
  }
  
  protocols <- list(
    distributed  = distr,
    concentrated = cnctr
  )
  
  if(input$strata=='overall'){
    
    results <- protocols %>% 
      map(
        .f = summarizer,
        data = abpm_long,
        input = input,
        truth = means_data[[input$data_lab]]
      ) %>% 
      bind_rows(.id='protocol')
    
  } else {
    
    results <- abpm_long %>% 
      split(f = .[[input$strata]]) %>% 
      map(
        ~protocols %>% 
          map(
            .f=summarizer,
            data = .x, 
            input = input, 
            truth = means_data[[input$data_lab]]
          ) %>% 
          bind_rows(.id='protocol')
      ) %>% 
      bind_rows(.id=input$strata)
    
  }
  
  saveRDS(
    results,
    paste0('Results/results_',paste(input,collapse='_'),'.RDS')
  )
  
  
}



# 
# oplt <-ggplot(results, aes(x=lvm,y=sbp,size=nht,label=lab,col=proto))+
#   geom_point()+
#   #geom_label_repel(show.legend = FALSE)+
#   theme_Publication()+
#   guides(size=FALSE,
#          colour = guide_legend(override.aes = list(size=10)))+
#   labs(size='',col='Sampling protocol',
#        y='Correlation with asleep \nsystolic blood pressure',
#        x='Correlation with left ventricular mass indexed to height')
# 
# if(input$strata!='none'){
#   
# }

#facet_wrap(~facet_var,scales='free')+

## Quality control
# View(analysis[analysis$subjid=='20407',])
# blk=analysis[analysis$subjid==analysis$subjid[1],]
# blk=analysis[analysis$subjid==analysis$subjid[1],]
# times=distr[[1]][[3]]
# impute=T
# dblk=data_blocks[[1]]
# proto=distr[[1]]
# t=proto[[1]]
















