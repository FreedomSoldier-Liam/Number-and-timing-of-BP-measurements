

tab_meansd <- function(variable,
                       data,
                       strat,
                       stratification_variable,
                       meta_variables,
                       include.pval=FALSE){
  
  mean_sd<-function(v){
    
    mu = mean(v, na.rm=TRUE)
    x_abs <- abs(mu)
    
    if(x_abs < 10) {
      dig = 2
    } else if(x_abs < 1000){
      dig = 1
    } else if(x_abs > 1000){
      dig = 0
    }
    
    mu %<>% round(dig) %>% format(nsmall=dig)
    sd = sd(v, na.rm=TRUE)%>%
      round(dig) %>% format(nsmall=dig)
    
    paste0(mu, ' (', sd,')')
    
  }
  
  r = gen_rowname.continuous_variable(meta_variables[[variable]])
  s = c(paste0("Overall<br/>(n = ", nrow(data),')'),strat$levels)
  t = tapply(data[[variable]], data[[strat$variable_name]],mean_sd)
  o = mean_sd(data[[variable]])
  
  if(include.pval){
    
    pval = if(length(strat$levels)==2){
      t.test(data[[variable]]~data[[strat$variable_name]]) %>% 
        use_series("p.value") %>%
        edit_pval()
    } else {
      as.formula(paste(variable,'~',strat$variable_name))%>%
        lm(data=data) %>% 
        anova() %>% 
        magrittr::extract(1,ncol(.)) %>% 
        edit_pval()
    }
    
    cbind(t(c(o,t)), p = pval) %>%
      set_rownames(r) %>%
      set_colnames(c(s, 'P Value'))
    
  } else {
    
    t(c(o,t))%>%
      set_rownames(r) %>%
      set_colnames(s)
    
  }
  
  
}
