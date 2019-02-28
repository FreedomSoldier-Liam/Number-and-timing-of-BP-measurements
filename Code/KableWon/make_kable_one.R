

make_kable_one <- function(tbl_one){
  
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
  
  return(kab)
  
}
