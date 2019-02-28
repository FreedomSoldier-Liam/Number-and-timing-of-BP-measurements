
longitudinal_analysis<-function(input,meta){
  
  exposure=input$exposure
  outcome<-input$outcome
  time=meta$root$variable[meta$root$group=='time']
  id=meta$root$variable[meta$root$group=='identification']
  
  all_times=levels(meta$data[[time]])
  
  model_vars=union(c(id,outcome,exposure),meta$root$variable[meta$root$model>0])
  model_frame=na.omit(meta$data[,model_vars])
  
  if(input$baseline.adj){
    
    bl_vals<-model_frame%>%plyr::ddply(.var=id,.fun=function(x){
      if(nrow(x)<2){
        return(NULL)
      } else {
        return(x[1,outcome])
      }
    })
    
    bl_mean<-mean(bl_vals$V1)
    
    model_frame%<>%plyr::ddply(.var=id,.fun = function(x){
      if(nrow(x)<2){
        return(NULL)
      } else {
        x$baseline_val=x[1,outcome]
        x[[outcome]]=x[2,outcome]
        return(x[-2,])
      }
    })
    meta$models%<>%map(.f=~gsub("~", "~ baseline_val +",.))
  }
  
  model_frame%<>%droplevels()
  
  nsubs=table(model_frame[[exposure]],model_frame[[time]])
  
  if(input$baseline.adj){
    cnames=colnames(nsubs)
    nsubs%<>%cbind(nsubs[,1])
    colnames(nsubs)=c('Baseline',cnames)
  }
  
  grps=meta$variables[[exposure]]$levels
  ngrps=length(grps)
  times = unique(model_frame[[time]])
  ntimes = length(times)
  nmdl=max(meta$root$model)
  
  names(model_frame)[names(model_frame)==input$outcome]='outcome'
  
  if(ntimes>1){
    
    models<-map(meta$models,.f=function(frm){
      lme4::lmer(as.formula(frm),data=model_frame)
    })
    
    nsubs%<>%
      as.data.frame%>%
      tidyr::spread(Var2,Freq)%>%
      mutate(Difference='--')
    names(nsubs)[1]=exposure
    
  } else {
    
    meta$models%<>%map(.f=~gsub(paste0('+ (1|',id,')'),"",.,fixed = TRUE))
    meta$models%<>%map(.f=~gsub(paste0('*',time),"",.,fixed = TRUE))
    models<-map(meta$models,.f=function(frm){
      lm(as.formula(frm),data=model_frame)
    })
    
    nsubs%<>%
      as.data.frame()%>%
      mutate(x=rownames(.))%>%
      mutate(Difference='--')%>%
      dplyr::select(x,everything())
    names(nsubs)[1]=exposure
    
  }
  
  require(r2glmm, quietly = TRUE)
  
  r2=map_chr(models,.f=function(m){
   format(round(r2beta(m,partial=F)$Rsq,2),nsmall=2)
  })
  
  smry=meta$data%>%dplyr::select_(.dots=list(id,time,outcome,exposure))%>%
    na.omit()%>%
    set_names(c("ID","t","o","e"))%>%
    dplyr::filter(ID%in%model_frame[[id]])%>%
    group_by(e,t)%>%
    dplyr::summarise(mn=mean(o),sd=sd(o),n=n())
  
  ybump=diff(range(meta$data[[outcome]],na.rm=TRUE))/5
  ymin=min(meta$data[[outcome]],na.rm=TRUE)-ybump
  
  descr_fig = meta$data[,c(id,time,outcome,exposure)]%>%
    droplevels%>%set_names(c("ID","t","o","e"))%>%
    dplyr::filter(ID%in%model_frame[[id]])%>%
    ggplot(aes(x=t,y=o))+
    geom_point(position=position_jitter(width = 1/15),col='grey60',size=0.9)+
    geom_pointrange(data=smry,aes(x=t,y=mn,ymin=mn-sd,ymax=mn+sd),
                    inherit.aes = FALSE, size=1,col='red',alpha=1)+
    geom_line(data=smry, aes(x=t,y=mn,group=1),col='red',size=1,
              alpha=0.50, linetype=1)+
    geom_text(data=smry,aes(x=t,y=ymin,label=paste("N =",n)),size=5)+
    facet_wrap(~e)+
    theme_Publication()+
    labs(x=paste0(meta$root$label[meta$root$variable==time],'\n'),
         y=paste0(meta$root$label[meta$root$variable==outcome],'\n'))
  
  diagn_figs<-map(models,.f=function(model){
    
    mf=model_frame%>%mutate(pred=predict(model))
    
    list(
      res.plot=mf%>%
        ggplot(aes(y=outcome-pred,x=pred,label=sid))+
        geom_point(alpha=0.33,col='grey50')+
        geom_hline(yintercept = 0)+
        geom_smooth(color='red',se=F)+
        geom_rug(sides='b')+
        labs(x='Predicted value',y='Residual value')+
        theme_Publication(),
      qq.plot=mf%>%
        ggplot(aes(sample=outcome-pred,label=sid))+
        stat_qq()+stat_qq_line()+
        theme_Publication()+
        labs(x='Quantiles of a standard normal distribution',
             y='Quantiles of residuals')
    )
    
  })
    
  kab = map(models,.f=function(model){
    
    if(inherits(model,'lmerMod')){
      
      null_form<-as.formula(paste(".~.-",exposure,":",time,"+(1|",id,")"))
      ovrl_form<-as.formula(paste("~",time))
      mean_form<-as.formula(paste("~",exposure,"|",time))
      intr_form<-as.formula(paste("~",exposure,"*",time))
      time_form<-as.formula(paste("~",time,"|",exposure))
      
      mnull=update(model,null_form,data=model@frame)
      intr_pval=edit_pval(pbkrtest::KRmodcomp(model,mnull)$test[1,'p.value'])
      
      
      #fitted_models<-list(model,model,model,mnull)
      
      EM = suppressMessages(
        list(mean_form,intr_form,time_form,ovrl_form)%>%
          map(.f=function(frm){
            
            mean_obj=emmeans::emmeans(object=model,frm)
            
            intr_obj=mean_obj%>%
              pairs(interaction=TRUE, reverse=T)%>%
              confint() %>% 
              as.data.frame() %>%
              mutate(tabval=paste0(adapt_round(estimate)," (",
                                   adapt_round(lower.CL),", ",
                                   adapt_round(upper.CL),")"))
            mean_obj%<>%summary()%>%as.data.frame()%>%
              mutate(tabval=paste0(adapt_round(emmean)," (",
                                   adapt_round(lower.CL),", ",
                                   adapt_round(upper.CL),")"))
            list(mean=mean_obj,intr=intr_obj)
          })%>%
          set_names(c("D1","D2","D3","D4"))
      )
      
      EM$D4$mean[,exposure]='Overall'
      
      tab1=rbind(EM$D1$mean[,c(exposure,time,'tabval')],
                 EM$D4$mean[,c(exposure,time,'tabval')])%>%
        tidyr::spread(time,tabval)
      
      if(input$baseline.adj){
        
        bl=expand.grid(
          c(grps,"Overall"),
          paste(adapt_round(bl_mean),"(ref)"))%>%
          set_names(c(exposure,'Baseline'))
        
        tab1=left_join(bl,tab1,by=exposure)%>%
          mutate_if(is.factor,as.character)%>%
          cbind(rep("",nrow(.)))%>%
          rbind(rep("",ncol(.)))
        
        tab1[nrow(tab1),1:(ncol(tab1)-1)]<-c("Difference",
                                             "0 (ref)",
                                             EM$D1$intr$tabval)
        
        tab1[,ncol(tab1)]=c(EM$D3$intr$tabval,
                            EM$D4$intr$tabval,
                            EM$D2$intr$tabval)
        
      } else {
        
        diffs_by_time=EM$D1$intr[,c(1:2,ncol(EM$D1$intr))]%>%
          tidyr::spread(exam,tabval)%>%
          set_names(names(tab1))
        
        tab1=tab1%>%
          rbind(diffs_by_time)%>%
          cbind(Difference=c(EM$D3$intr$tabval,
                             EM$D4$intr$tabval,
                             EM$D2$intr$tabval))
        
        # times<-meta$variables[[time]]$levels
        # bl<-times[1]; fu<-times[length(times)]
        # kp<-paste(fu,'-',bl)
        # 
        # tab1[,ncol(tab1)]=c(EM$D3$intr$tabval[EM$D3$intr[,1]==kp],
        #                     EM$D4$intr$tabval[EM$D4$intr[,1]==kp],
        #                     EM$D2$intr$tabval[EM$D2$intr[,2]==kp])
        
      }
      
      #names(tab1)[ncol(tab1)]='Difference'
      
      tab1[[exposure]]%<>%factor(levels=c(
        "Overall",levels(model_frame[[exposure]]),
        levels(EM$D2$intr[,1])))
      
      tab1%<>%dplyr::arrange_(.dots=exposure)
      
      tab1
      
    } else if(inherits(model,'lm')){
      
      # tst=anova(model)
      # intr_pval = tst[rownames(tst)==exposure,ncol(tst)]%>%edit_pval()
      
      bl=expand.grid(
        c(grps,"Overall"),
        paste(adapt_round(bl_mean),"(ref)"))%>%
        set_names(c(exposure,'Baseline'))
      
      mean_obj=as.formula(paste("~",exposure))%>%
        emmeans::emmeans(object=model, frm)
      
      ovrl_obj=as.formula(paste("~1"))%>%
        emmeans::emmeans(object=model, frm)%>%
        confint()%>%as.data.frame()
      
      names(ovrl_obj)[1]=exposure
      ovrl_obj[,exposure]="Overall"
      
      EM=list(
        intr=mean_obj%>%
          pairs(interaction=TRUE, reverse=TRUE)%>%
          confint() %>% 
          as.data.frame() %>%
          mutate(tabval=paste0(adapt_round(estimate)," (",
                               adapt_round(lower.CL),", ",
                               adapt_round(upper.CL),")")),
        mean=mean_obj%>%
          summary()%>%as.data.frame()%>%
          mutate(`Follow up`=paste0(adapt_round(emmean)," (",
                                    adapt_round(lower.CL),", ",
                                    adapt_round(upper.CL),")"),
                 Difference=paste0(adapt_round(emmean-bl_mean)," (",
                                   adapt_round(lower.CL-bl_mean),", ",
                                   adapt_round(upper.CL-bl_mean),")")),
        ovrl=ovrl_obj%>%
          mutate(`Follow up`=paste0(adapt_round(emmean)," (",
                                    adapt_round(lower.CL),", ",
                                    adapt_round(upper.CL),")"),
                 Difference=paste0(adapt_round(emmean-bl_mean)," (",
                                   adapt_round(lower.CL-bl_mean),", ",
                                   adapt_round(upper.CL-bl_mean),")"))
      )
      
      merge_vars=c(exposure,'Follow up','Difference')
      
      tab1=bl %>%
        merge(rbind(EM$mean,EM$ovrl)[,merge_vars],by=exposure)%>%
        mutate_if(is.factor,as.character)
      
      for(i in 1:nrow(EM$intr)){
        tab1%<>%rbind(rep("",ncol(.)))
        tab1[nrow(tab1),]=c(paste(EM$intr[i,1]),
                            "0 (ref)",
                            rep(EM$intr[i,'tabval'],2))
      }
      
      tab1[[exposure]]%<>%factor(levels=c(
        "Overall",levels(model_frame[[exposure]]),levels(EM$intr[[1]])))
      
      names(tab1)[names(tab1)=='Follow up']=paste(all_times[-1])
      
      tab1%<>%dplyr::arrange_(.dots=exposure)
      
      tab1
    }
    
  }) %>% reduce(rbind)
  
  names(nsubs)=names(kab)
  kab=rbind(nsubs,kab)
  
  header=c("",ncol(kab)-1)
  names(header)=c("",meta$root$label[meta$root$variable==input$outcome])
  names(kab)[1]=''
  
  kab %<>%
    knitr::kable(
      format='html',align=rep('c', ncol(.)), 
      booktabs = TRUE,escape = FALSE)%>%
    kable_styling(bootstrap_options=c("striped","hover"),full_width=T)%>% 
    column_spec(1, width = "12cm")%>%
    add_header_above(header)%>%
    group_rows(group_label='Group Size',start_row=1,end_row=ngrps)
  
  ncmp=choose(ngrps,2)
  
  breaks=ngrps+(0:nmdl)*(1+ncmp+ngrps)
  group_labels=paste0("Model ",1:nmdl,"[note]")
  
  for(i in 1:nmdl){
    kab %<>% group_rows(group_label = group_labels[i],
                        start_row = breaks[i]+1,
                        end_row = breaks[i+1])
  }
  
  notes <- map(1:nmdl,.f=function(i){
    paste0("Model ",i,' included ', 
          ifelse(i==1,'adjustment','additional adjustment'),' for ',
          list_elements(tolower(meta$root$label[meta$root$model==i])),
          " (R-squared =", r2[i],')')
  })
  
  kab%<>%add_footnote(notes)
  
  list(kab=kab,descr_fig=descr_fig,diagn_figs=diagn_figs)
  
}