
make_table_one <- function(meta, time=NULL,
                           include.pval=TRUE, 
                           include.freq=FALSE,
                           caption=NULL,
                           as.kable=FALSE){
  
  stratification_variable=meta$exposure
  
  is.longitudinal_data<-!is.null(time)
  
  if(is.longitudinal_data){
    time_variable=meta$root$variable[meta$root$group=='time']
  }

  id=meta$root$variable[meta$root$group=='identification']
  meta_variables=meta$variables
  
  data=meta$data
  data[[id]]=NULL
  
  if(is.longitudinal_data){
    data=data[data[[time_variable]]==time,]
  }

  
  # Remove meta variables that aren't in the given data set
  for(m in names(meta_variables)){
    if(!(m%in%names(data))) meta_variables[[m]]=NULL
    if("group"%in%names(meta_variables[[m]])){
      if(meta_variables[[m]]$group=='time') meta_variables[[m]]=NULL
    }
  }
  
  # Pull the stratifying variable out of the meta variable list
  strat = meta_variables[[stratification_variable]]
  strat$variable_name=stratification_variable
  
  # Create header labels using stratification variable
  t = table(data[[strat$variable_name]])
  strat$levels %<>% paste0('<br/>(n = ', t, ')')
  
  # Add a footnote if there is one printed in the meta file
  if(!is.null(strat$footnote)){
    if(!is.na(strat$footnote)){
      if(strat$footnote!=''){
        strat$label%<>%paste0("[note]")
      } 
    }
  }

  # Identify variable types in the data
  var_class=unique(map_chr(meta_variables,~class(.)))
  var_types=map(var_class,.f=function(vc){
    map_lgl(meta_variables,~class(.)==vc)%>%which()%>%names()
  }) %>% set_names(var_class)
  
  # Construct footnote for bottom of the table
  n_str <- if (include.freq) 'count (%)' else 'percentage'
  m_str <- 'mean (SD)'
  
  rmv_blank=function(x){
    
    if(any(!is.na(x))){
      if(any(x=='NA')){
        x=x[-which(x%in%c('NA'))]
      }
      if(any(x=='') | any(x==" ")){
        x=x[-which(x%in%c('',' '))]
      }
    }
    x
  } 
  
  abbrv=map_chr(meta_variables,.f=~.$abbreviation)%>%
    unique()%>%rmv_blank()
  
  if(length(abbrv)>0){
    abbrv%<>%str_split(';')%>%
      unlist()%>%trimws()%>%sort()%>%paste(collapse='; ')
  } else {
    abbrv=''
  }
  
  c1=length(var_types$categorical_variable)>0
  c2=length(var_types$continuous_variable)>0
  
  if(c1&c2){
    notes=paste0('Table values are ', m_str, ' or ', n_str,'.')
  } else if(c1&!c2){
    notes = paste0('Table values are ', n_str,'.')
  } else if(!c1&c2){
    notes = paste0('Table values are ', m_str,'.')
  }
  
  # add footnotes to the notes if user specified notes in meta data
  notes = c(notes,map_chr(meta_variables, ~ paste0(.$footnote))) %>%
    c(abbrv)%>%
    rmv_blank()
  
  omit_indx=grep(stratification_variable,names(meta_variables))
  
  output <- meta_variables[-omit_indx] %>% 
    map2(names(meta_variables)[-omit_indx],.f=function(v,l){
      
      if(inherits(v,'continuous_variable')){
        tab_meansd(variable=l,
                   data=data,
                   strat=strat,
                   stratification_variable=stratification_variable,
                   meta_variables=meta_variables,
                   include.pval=include.pval)
      } else if (inherits(v,'categorical_variable')){
        tab_colperc(variable=l,
                    data=data,
                    strat=strat,
                    stratification_variable=stratification_variable,
                    meta_variables=meta_variables,
                    include.freq=include.freq,
                    include.pval=include.pval)
      } 
      
    }) %>% reduce(rbind)
  
  if(is.null(caption)){
    caption <- paste0(
      'Participant characteristics[note] overall and stratified by ', 
      tolower(strat$label),'.')
  }
  
  tbl_one = structure(
    list(
      table = output,
      meta_variables=meta_variables,
      strat=strat,
      include.pval=include.pval,
      caption=caption,
      notes=notes,
      abbrv=abbrv,
      class = 'Table_One'
    ))
  
  if(as.kable){
    omit_indx=grep(tbl_one$strat$variable_name,names(tbl_one$meta_variables))
    groups = map_chr(tbl_one$meta_variables[-omit_indx],~paste(.$group))
    group_labels=unique(groups)
    breaks = c(0)
    nlevels=map_dbl(tbl_one$meta_variables[-omit_indx],~length(.$levels))
    nlevels[nlevels<3]=0
    
    for(i in 1:length(group_labels)){
      breaks[i+1] = breaks[i] + sum(
        map_dbl(tbl_one$meta_variable[-omit_indx],
                ~.$group==group_labels[i])*(1+nlevels)
      )
    }
    
    kab = suppressWarnings(
      knitr::kable(
        tbl_one$table,format='html',align=rep('c', ncol(tbl_one$table)), 
        caption = tbl_one$caption, booktabs = TRUE,escape = FALSE)%>%
        kable_styling(bootstrap_options=c("striped","hover"),full_width=T)%>% 
        column_spec(1, width = "6cm")%>% 
        add_indent(which(grepl("    ",rownames(tbl_one$table))))%>% 
        add_footnote(tbl_one$notes,notation='symbol',threeparttable=T)
    )
    
    for(i in 1:length(group_labels)){
      kab %<>% group_rows(group_label = group_labels[i],
                          start_row = breaks[i]+1,
                          end_row = breaks[i+1])
    }
    
    tbl_one=kab
    
  }
  
  return(tbl_one)
  
}
