rm(list=ls())
setwd("~/Dropbox/MOST_NSC/107年計畫107年八月/snp/proposal/plot/")
#install.packages("PIGShift")
library(PIGShift)
library(ape)
library(phangorn)
library(TreeSim)
#library(xtable)
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV_test.r")
#load("plotcluster.RData")
size<-6
tree <-  sim.bd.taxa.age(n=size, numbsim=1, lambda=2, mu=1, frac = 0.5, age=1, mrca = TRUE)[[1]]
#tree2 <- DV(k=2,tree=tree)$V
#tree3 <- DV(k=3,tree=tree)$V
#jpeg("alp_cluster.jpeg")
par(mfrow=c(1,3))
tree0<-vcv2phylo(vcv(tree))
tree0$tip.label<-LETTERS[1:size]
plot(tree0, main=expression(paste(alpha, " = 0",sep="")), edge.width=3,cex.main=3,cex=2, direction="upward")

theta<-2
V2<-OU.vcv(tree,theta=theta)
#tree2<-vcv2phylo(round(VaF(alpha=2,V=V),2))
tree2<-vcv2phylo(V2)
tree2$tip.label<-LETTERS[1:size]
plot(tree2, main=expression(paste(alpha, " = ",1, sep="")), edge.width=3,cex.main=3,cex=2, direction="upward")

theta<-10
V3<-OU.vcv(tree,theta=theta)
#tree2<-vcv2phylo(round(VaF(alpha=2,V=V),2))
tree3<-vcv2phylo(V3)
tree3$tip.label<-LETTERS[1:size]
plot(tree3, main=expression(paste(alpha, " = ",10, sep="")), edge.width=3,cex.main=3,cex=2,direction="upward")
