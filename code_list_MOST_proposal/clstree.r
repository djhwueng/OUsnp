#https://stackoverflow.com/questions/33730729/how-display-length-of-branches-in-phylogenetic-tree

rm(list=ls())
setwd("~/Dropbox/MOST_NSC/107年計畫107年八月/snp/proposal/Latex/")
library(ape)
library(phangorn)
library(TreeSim)
library(xtable)
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV_test.r")
size<-6
tree <-  sim.bd.taxa.age(n=size, numbsim=1, lambda=1, mu=1, frac = 0.5, age=1, mrca = TRUE)[[1]]
tree2 <- DV(k=2,tree=tree)$V
tree3 <- DV(k=3,tree=tree)$V
#plot(tree)
#load("plotcluster.RData")
#plot(tree)

#jpeg("clstree.jpeg")
par(mfrow=c(1,3))
ck0tree <- vcv2phylo(vcv(tree))
ck0tree$tip.label<-LETTERS[1:size]
plot(ck0tree,edge.width=3,main=expression(paste(k, " = ", 0,sep="")),cex.main=4,cex=2,direction="upward")

ck2tree <- vcv2phylo(tree2)
ck2tree$tip.label<-LETTERS[1:size]
plot(ck2tree,edge.width=3,main=expression(paste(k, " = ", 2,sep="")),cex.main=4,cex=2,direction="upward")

ck3tree <- vcv2phylo(tree3)
ck3tree$tip.label<-LETTERS[1:size]
plot(ck3tree,edge.width=3,main=expression(paste(k, " = ", 3,sep="")),cex.main=4,cex=2,direction="upward")
#dev.off()

save.image("plotcluster.RData")

# #https://tex.stackexchange.com/questions/55418/how-to-print-a-latex-matrix-using-xtable-in-r
 xtable(round(vcv(ck0tree),2))
 xtable(round(vcv(ck2tree),2))
 xtable(round(vcv(ck3tree),2))
#
#  a<-matrix(rnorm(25),5,5)
#  x<-xtable(a,align=rep("",ncol(a)+1))
#  print(x,floating=FALSE,tabular.environment="bmatrix",hline.after=NULL,include.rownames=FALSE,include.colnames=FALSE)
# #
#  A<-matrix(NA,nrow=1,ncol=3)
#  xtable(A)
