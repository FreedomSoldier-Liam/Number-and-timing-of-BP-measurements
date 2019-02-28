
create_meta_object <- function(real_data,meta_data,exposure){
  
  meta_data %<>% dplyr::filter(variable%in%names(real_data))
  real_data=real_data[,paste(meta_data$variable)]
  
  meta_variables = vector(mode='list',length=ncol(real_data)) %>% 
    set_names(names(real_data))
  
  for(i in 1:nrow(meta_data)){
    
    mta=meta_data[i,]
    
    ## Variable processing -- CATEGORICAL
    if(mta['type']=='factor'){
      
      key = mta$key %>% 
        stringr::str_split(pattern=',') %>%
        unlist() %>%
        trimws() %>%
        map(.f=~unlist(str_split(.,pattern=' = ')))
      
      keyvals <- map_chr(key,~.[2])
      names(keyvals) <- map_chr(key,~.[1])
      
      real_data[,paste(mta$variable)] %<>% 
        factor() %>% 
        plyr::revalue(replace=keyvals)
      
      meta_variables[[i]] = structure(
        list(label = mta$label, 
             levels = keyvals%>%set_names(NULL),
             footnote=mta$note,
             abbreviations=mta$abbreviations,
             group=ifelse(i==exposure,'main exposure',mta$group),
             model=mta$model),
        class = 'categorical_variable'
      )
      
    }
    
    ## Variable processing -- CONTINUOUS
    if(mta$type%in%c('integer','numeric')){
      
      meta_variables[[i]] = structure(
        list(label= mta$label, 
             unit = mta$key, 
             footnote=mta$note,
             abbreviations=mta$abbreviations,
             group = mta$group,
             model = mta$model),
        class = 'continuous_variable'
      )
    }
    
  }
  
  model_formulas=vector(mode='list',length=max(meta_data$model))
  
  if('time'%in%meta_data$group){
    time=meta_data$variable[meta_data$group=='time']
    model_formulas[[1]]=paste0('outcome ~ ',exposure,'*',time)
  } else {
    model_formulas[[1]]=paste0('outcome ~ ', exposure)
  }
  
  for(i in 2:length(model_formulas)){
    model_formulas[[i]]=paste(model_formulas[[i-1]],'+',
                              paste(meta_data$variable[meta_data$model==i],
                                    collapse=' + '))
  }
  
  if('time'%in%meta_data$group){
    model_formulas%<>%map(~paste0(
      .," + (1|",meta_data$variable[meta_data$group=='identification'],")"))
  }
  
  
  list(root=meta_data,
       data=real_data,
       variables=meta_variables,
       exposure=exposure,
       models=model_formulas)
  
}
