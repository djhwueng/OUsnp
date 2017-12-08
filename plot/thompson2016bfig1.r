#http://www.phytools.org/eqg2015/asr.html
#http://schmitzlab.info/phylo.html
rm(list=ls())
library(ape)
mini.phy <- read.tree(text = "((((A,B), C), (D,E)),F);") #defining a tree in bracket form
#mini.phy$edge
edge.co <- c("red","black","black","black","black","black","red","red","red","black" )
tipcol<-c("black", "white", "white", "red", "red", "black")

setwd("~/Dropbox/MOST_NSC/107年計畫107年八月/snp/proposal/Latex/")
png("T2014Fig.png")
plot(mini.phy,edge.col = edge.co, adj=0, label.offset=0.75,edge.width=5)
tiplabels(pch=c(22,2,2,21,21,22), col="black", adj=1, bg=tipcol, cex=c(1,1,2,3,3,5))
axisPhylo(backward=F)
dev.off()
?tiplabels
