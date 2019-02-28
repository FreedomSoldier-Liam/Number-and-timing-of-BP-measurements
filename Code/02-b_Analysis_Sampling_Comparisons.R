rm(list=ls())

theme_Publication <- function(base_size=16) {
  
  require(ggthemes)
  
  (theme_foundation(base_size=base_size)+ 
      theme(plot.title = element_text(face = "bold",
                                      size = rel(1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.border = element_rect(colour = 'black'),
            legend.key.size = unit(3,"line"),
            legend.key.height = unit(3,"line"),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size=rel(1)),
            axis.line = element_blank(), # element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(), #element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title = element_text(face="italic"),
            legend.text = element_text(size=rel(1)),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="black",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
      ))
  
}
mae<-function(x,y) mean(abs(x-y),na.rm=TRUE)

library(tidyverse)
library(magrittr)
library(fmsb)
library(BlandAltmanLeh)

## Distributed protocol times

lwr_time=1
upr_time=5
ntimes=list(2,3,4)

distr=map(ntimes, ~combn(lwr_time:upr_time,.))%>%
  map(.f=function(mat){map(1:ncol(mat),~mat[,.])})%>%
  set_names(paste("msr",ntimes,'times',sep='_'))

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
  data_lab=c("CARDIA","JHS"),
  impute=c(FALSE,TRUE),
  strata=c("overall",'sex','race'),
  tvarbl=c('tss','tsm'),
  stringsAsFactors = FALSE
)

take_out=which(inputs$strata=='race'&inputs$data_lab=='JHS')
inputs=inputs[-take_out,]

i=1

for(i in 1:nrow(inputs)){
  
  input=inputs[i,]
  
  analysis=readRDS(paste0("Datasets/",input$data_lab,"_abpm.RDS"))
  analysis$time=analysis[[input$tvarbl]]
  
  if(input$impute){
    analysis%<>%dplyr::mutate(sbp=sbp_imp,dbp=dbp_imp)
  }
  
  # blk=analysis[analysis$subjid==analysis$subjid[1],]
  # times=distr[[1]][[2]]
  # impute=T
  
  sampler<-function(blk,times,impute=F){
    
    indx<-map_int(times,~which.min(abs(blk$time-.)))
    
    if(any(blk$time[indx]-times>1.5)){
      print(blk$subjid[1])
      return(NULL)
    } else {
      blk[indx,]
    }
    
  }
  
  # dblk=analysis
  # proto=cnctr$msr_3_times
  
  summarizer<-function(proto,dblk){
    
    map(proto,.f=function(t){
      
      smp=dblk%>%
        plyr::ddply('subjid',sampler,times=t,input$impute)%>%
        group_by(subjid)%>%
        dplyr::summarise(slp_sbp_est=mean(sbp),
                         slp_dbp_est=mean(dbp),
                         slp_sbp=slp_sbp[1],
                         slp_dbp=slp_dbp[1],
                         nht=nht[1])%>%
        mutate(nht2=ifelse(slp_sbp_est>=120|slp_dbp_est>=70,1,0))
      
      tlab=paste(t,collapse='_')
      tlab=gsub('.','deci',tlab,fixed = TRUE)
      fname=paste0(paste(input,collapse='_'),"_",tlab,".RDS")
      saveRDS(smp,file.path('Results','Full Results',fname))
      
      # p.sbp=bland.altman.plot(smp$slp_sbp_est,smp$slp_sbp,
      #                         graph.sys="ggplot2")+
      #   theme_Publication()+
      #   labs(x='Mean of measurements',
      #        y="Difference in measurements",
      #        title='Systolic blood pressure')
      # 
      # p.dbp=bland.altman.plot(smp$slp_dbp_est,smp$slp_dbp,
      #                         graph.sys="ggplot2")+
      #   theme_Publication()+
      #   labs(x='Mean of measurements',
      #        y="Difference in measurements",
      #        title='Diastolic blood pressure')
      # 
      # p.all=cowplot::plot_grid(p.sbp,p.dbp,nrow=1)
      
      data.frame(
        sbp=mae(smp$slp_sbp_est,smp$slp_sbp),
        dbp=mae(smp$slp_dbp_est,smp$slp_dbp),
        #lvm=cor(smp$slp_sbp_est,smp$lvm),
        agr=mean(smp$nht==smp$nht2),
        nht=Kappa.test(x=smp$nht,y=smp$nht2)$Result$estimate,
        lab=paste(t,collapse='-')
      )
      
    }) %>% 
      reduce(rbind)
    
  }
  
  polisher<-function(dblk){
    
    list(
      
      distr=map(distr, .f=summarizer,dblk)%>%
        reduce(rbind),
      
      cnctr=map(list(cnctr[[input$data_lab]]),.f=summarizer,dblk)
      
    ) %>%
      map2(names(.),.f=function(res,proto){
        res%>%data.frame()%>%
          mutate(lab=as.character(lab),proto=proto)})%>%
      reduce(rbind)
    
  }
  
  if(input$strata=='overall'){
    results=polisher(analysis)
  } else {
    data_blocks<-analysis%>%split(f=.[[input$strata]])
    results=map(data_blocks,.f=polisher)%>% 
      map2(names(.),.f=function(dblk,dlab){
        dblk%>%mutate(facet_var=dlab)
      })%>%
      reduce(rbind)
  }
  
  results %>% 
    dplyr::mutate(proto=factor(proto,levels=c("cnctr","distr"),
                               labels=c("Concentrated","Distributed")))%>%
    saveRDS(paste0('Results/results_',paste(input,collapse='_'),'.RDS'))
  
  
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
















