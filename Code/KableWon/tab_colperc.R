
tab_colperc <- function(variable,
                        data,
                        strat,
                        stratification_variable,
                        meta_variables,
                        include.freq=FALSE,
                        include.pval=FALSE){
  
  indent <- function(x) paste('   ', x)
  
  m = meta_variables[[variable]]
  t = table(data[[variable]], data[[stratification_variable]])
  o = table(data[[variable]])
  pt = format(round(100*prop.table(t, margin = 2),1),nsmall=1) 
  po = format(round(100*prop.table(o),1),nsmall=1)
  
  t=cbind(o,t)
  p=cbind(po,pt)
  
  r = gen_rowname.categorical_variable(meta_variables[[variable]])
  s = c(paste0("Overall<br/>(n = ", nrow(data),')'),strat$levels)
  
  ## Categorical variable has 2 levels
  if(nrow(t)==2){
    
    t2=t[2,] %>% trimws(); p2=p[2,] %>% trimws()
    
    if(include.freq){
      vals <- paste0(t2, ' (',p2, '%)')
    } else {
      vals <- paste0(p2,'%')
    }
    
    res = matrix(vals, ncol = length(t2)) %>% 
      set_colnames(s) %>% set_rownames(r)
    
    if(include.pval){
      chi_pval = chisq.test(t)$p.value %>% edit_pval()
      chi_pvec = c(chi_pval)
    }
    
  } 
  
  ## Categorical variable has >2 levels
  if(nrow(t)>2){
    
    if(include.freq){
      vals <- paste0(t, ' (', p, '%)')
    } else {
      vals <- paste0(p,'%')
    }
    
    res = c("") %>% 
      rbind(matrix(vals, ncol = ncol(t))) %>%
      set_colnames(s) %>%
      set_rownames(c(r,indent(paste0(
        meta_variables[[variable]]$levels
      ))))
    
    if(include.pval){
      chi_pval = chisq.test(t)$p.value %>% edit_pval()
      chi_pvec = c(chi_pval, rep("",nrow(t)))
    }
    
  }
  
  if(include.pval){
    cbind(res, 'P Value' = chi_pvec)
  } else {
    res
  }
  
}