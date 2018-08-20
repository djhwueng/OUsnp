# http://www.gwaspi.org/?page_id=671
# To install core packages:
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# To install specific packages
# biocLite("snpStats")
# Read a PLINK binary data file as a SnpMatrix
# Description
# The package PLINK saves genome-wide association data in groups of three files, with the extensions .bed, .bim, and .fam. This function reads these files and creates an object of class "SnpMatrix"
rm(list=ls())
#setwd("~/Documents/drosophila")
setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/drosophila")
# install.packages("ape")
# install.packages("phytools")
# install.packages("geiger")
# install.packages("phyclust")
# install.packages("corpcor")
# install.packages("phylobase")
# install.packages("gap")
# install.packages("TreeSim")
# install.packages("OUwie")

library(snpStats)
library(nlme)
library(geiger)
library(ape)
library(phytools)
#source("~/Dropbox/YiWeiHsu/Rcode/main.code/DV_test.r")
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV_test.r")
#source("~/Dropbox/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
# snpmtx<-read.plink("dgrp2.bed","dgrp2.bim.txt","dgrp2.fam.txt", na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = NULL)
# snpgeno<-as.data.frame(snpmtx$genotypes)
# snpmap<-snpmtx$map
# rm(snpmtx)

load("snpdata.RData")
dim(snpgeno) #  4 million genome x 205  flies 

chr6g<-snpgeno[,4427848:4438427]
chr6m<-snpmap[snpmap$chromosome==6,]
chr6gp<-matrix(as.numeric(t(chr6g)),ncol=205)
chr6gp<-t(chr6gp)
rownames(chr6gp)<-rownames(chr6g)
colnames(chr6gp)<-chr6m$position
table(chr6gp)
chr6gp[chr6gp==0]<-1
chr6gp[chr6gp==3]<-0

write.nexus.data(chr6gp,file="drs6.nex")
system("sed 's/DNA/standard/g' drs6.nex >drs6.1.nex")
system("paup -n get.tree.snp.drs6.nex >/dev/null")
system("grep 'PAUP_1' drs6.tre | sed 's/^.*: //' > drs6.1.tre ")
temp.tree<-read.table(file="drs6.1.tre")
tree<-unlist(temp.tree[5])
tree<-as.character(tree)
write(tree,file="treedrs6.txt")
tree<-read.newick(file="treedrs6.txt")
writeNexus(tree,file="drs6.nex")


drs6<-read.nexus("drs6.nex")
plot(drs6)
#trait<-read.csv("/Users/yiweihsu/Documents/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)
trait<-read.csv("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)

head(trait)
rownames(trait)<-trait$V1
head(trait)
obj<-name.check(drs6,trait)
obj



dim(trait)
drs6.trait<-drop.tip(drs6,obj$tree_not_data)
name.check(drs6.trait,trait)

del1<-which(rownames(chr6gp)=="line_256")
del2<-which(rownames(chr6gp)=="line_509")
delset<-c(del1,del2)

chr6gprm<-chr6gp[-delset,]

bm<-corBrownian(1,drs6.trait)
p.value.array_6th<-array(0,c(dim(chr6gprm)[2]))
for(snpIndex in 1:dim(chr6gprm)[2]){
  print(snpIndex)
  YX<-cbind(trait$V2,chr6gprm[,snpIndex])
  colnames(YX)<-c("Y","X")
  YX<-as.data.frame(YX)
  try(model01<-gls(Y~X,data=YX, correlation=bm))
  if(exists("model01")){
  summarymodel<-summary(model01)
  p.value.array_6th[snpIndex]<-summarymodel$tTable[,4][2]
  }else{p.value.array_6th[snpIndex]<-NA}
  rm(model01)
  }
plot(-log10(p.value.array_6th))
save.image("drs6th.RData")

jpeg("drs6th.jpeg")
plot(-log10(p.value.array_6th))
dev.off()
