rm(list=ls())
library(phytools)
setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/micezhangdataset/")
load("realtree.RData")

y5<-list(phy3=phy3, phy4=phy4)
class(y5) <- "multiPhylo"
contree5<-ls.consensus(y5,tol=1e-3)
save.image("contree5.RData")
