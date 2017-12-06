rm(list=ls())
library(phytools)
setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/micezhangdataset/")
load("realtree.RData")

y1<-list(phy1=phy1, phy2=phy2)
class(y1) <- "multiPhylo"
contree1<-ls.consensus(y1,tol=1e-3)
save.image("contree1.RData")
