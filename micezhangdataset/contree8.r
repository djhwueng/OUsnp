rm(list=ls())
library(phytools)
setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/micezhangdataset/")
load("realtree.RData")

y8<-list(phy5=phy5, phy6=phy6)
class(y8) <- "multiPhylo"
contree8<-ls.consensus(y8,tol=1e-3)
save.image("contree8.RData")
