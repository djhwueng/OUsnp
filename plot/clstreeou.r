rm(list=ls())
setwd("~/Dropbox/MOST_NSC/107年計畫107年八月/snp/proposal/plot/")
library(ape)
library(phangorn)
library(TreeSim)
library(xtable)
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV_test.r")
#load("plotcluster.RData")
size<-6
tree <-  sim.bd.taxa.age(n=size, numbsim=1, lambda=2, mu=1, frac = 0.5, age=1, mrca = TRUE)[[1]]
#tree2 <- DV(k=2,tree=tree)$V
#tree3 <- DV(k=3,tree=tree)$V
V<-vcv(tree)
VaF<-function(alpha,V=V){
  A<-exp(-2*alpha*(1-V))
  B<- (1-exp(-2*alpha*V))/(2*alpha)
  return(A*B)
  }

jpeg("alp_cluster.jpeg")
par(mfrow=c(1,2))
tree0<-vcv2phylo(vcv(tree))
tree0$tip.label<-LETTERS[1:size]
plot(tree0, main=expression(paste(alpha, " = 0",sep="")), edge.width=5,cex.main=3,cex=3 )

tree2<-vcv2phylo(round(VaF(alpha=2,V=V),2))
tree2$tip.label<-LETTERS[1:size]
plot(tree2, main=expression(paste(alpha, " = 2",sep="")), edge.width=5,cex.main=3,cex=3 )
dev.off()
