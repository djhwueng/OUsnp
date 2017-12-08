#http://schmitzlab.info/visual1.html
rm(list=ls())
library(ape)
library(geiger)
library(phytools)
#library(treeman)
#install.packages("treeman")
tree<-read.tree(url("http://schmitzlab.info/visual.tree.tre"))
tree
visual.data<-read.csv(url("http://schmitzlab.info/visual.data.csv"),header=TRUE)
head(visual.data)
rownames(visual.data)<-visual.data$taxon


tip.to.drop<-paste("tip",which(visual.data$group=="unknown"),sep="")
tree1<-drop.tip(tree,tip.to.drop)
visual.data1<-visual.data[-which(visual.data$group=="unknown"),]
compare<-treedata(tree1,visual.data1,sort=TRUE)

tree<-compare$phy
visual.data<-as.data.frame(compare$data)
length(tree$tip.label)

z<-as.factor(visual.data$group); names(z)<-tree$tip.label

mycol<-character(length(z))
mycol[visual.data$groups=="cathemeral"]<-"blue"
mycol[visual.data$groups=="diurnal"]<-"yellow"
mycol[visual.data$groups=="nocturnal"]<-"red"
#mycol[visual.data$groups=="unknown"]<-"green"
par(mar=c(2,0,0,1))
setwd("~/Dropbox/MOST_NSC/107年計畫107年八月/snp/proposal/plot/")
jpeg("treesnp.jpeg")
plot(tree, show.tip.label=FALSE)#,direction="upwards")
tiplabels(pch=21,col="black",bg=mycol,cex=1.5)
dev.off()
#nodelabels()
