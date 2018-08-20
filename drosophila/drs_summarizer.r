rm(list=ls())
setwd("/Users/yiweihsu/Documents/drosophila")
rdata<-c("drs1th.RData","drs2th.RData","drs3th.RData","drs4th.RData","drs5th.RData","drs6th.RData")
par(mfrow=c(2,3))
for(dIndex in 1:length(rdata)){
  #dIndex<-1
  load(rdata[dIndex])
  p_value<-get(paste("p.value.array_",dIndex,"th",sep=""))
  plot(-log10(p_value),ylab= "-log10 p", xlab="Drosophila position",main= paste("chromosome ", dIndex, sep=""))
  }