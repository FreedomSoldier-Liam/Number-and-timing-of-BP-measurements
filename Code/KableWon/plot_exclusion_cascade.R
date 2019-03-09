
library(diagram)
library(tidyverse)


plot_exclusion_cascade <- function(nobs, labs, txt_expand){
  
  par(mar = c(1, 1, 1, 1))
  nrows=length(nobs)
  openplotmat()
  pos<-coordinates(rep(1,nrows))
  
  for(i in 1:nrows){
    
    if(i<nrows){
      Arrows(x0=pos[i,1],y0=pos[i,2]-1/(nrows*3.8),
             x1=pos[i+1,1],y1=pos[i+1,2]+1/(nrows*3.8),lwd=2)
    }
    textrect(mid = pos[i,], radx = 1/2, rady = 1/(nrows*4),
             lab = paste0(labs[[i]],'\n N = ',format(nobs[[i]],big.mark=',')), 
             box.col='white', col='black',
             cex = txt_expand, shadow.size = 0.001)
    
  }
}




