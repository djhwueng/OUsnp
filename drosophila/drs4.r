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

chr4g<-snpgeno[,2840330:3827021]
chr4m<-snpmap[snpmap$chromosome==4,]
chr4gp<-matrix(as.numeric(t(chr4g)),ncol=205)
chr4gp<-t(chr4gp)
rownames(chr4gp)<-rownames(chr4g)
colnames(chr4gp)<-chr4m$position
table(chr4gp)
chr4gp[chr4gp==0]<-1
chr4gp[chr4gp==3]<-0

write.nexus.data(chr4gp,file="drs4.nex")
system("sed 's/DNA/standard/g' drs4.nex >drs4.1.nex")
system("paup -n get.tree.snp.drs4.nex >/dev/null")
system("grep 'PAUP_1' drs4.tre | sed 's/^.*: //' > drs4.1.tre ")
temp.tree<-read.table(file="drs4.1.tre")
tree<-unlist(temp.tree[5])
tree<-as.character(tree)
write(tree,file="treedrs4.txt")
tree<-read.newick(file="treedrs4.txt")
writeNexus(tree,file="drs4.nex")


drs4<-read.nexus("drs4.nex")
plot(drs4)
#trait<-read.csv("/Users/yiweihsu/Documents/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)
trait<-read.csv("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)
head(trait)
rownames(trait)<-trait$V1
head(trait)
obj<-name.check(drs4,trait)
obj
dim(trait)
drs4.trait<-drop.tip(drs4,obj$tree_not_data)
name.check(drs4.trait,trait)


del1<-which(rownames(chr4gp)=="line_256")
del2<-which(rownames(chr4gp)=="line_509")
delset<-c(del1,del2)

chr4gprm<-chr4gp[-delset,]
bm<-corBrownian(1,drs4.trait)
p.value.array_4th<-array(0,c(dim(chr4gprm)[2]))
for(snpIndex in 1:dim(chr4gprm)[2]){
  print(snpIndex)
  YX<-cbind(trait$V2,chr4gprm[,snpIndex])
  colnames(YX)<-c("Y","X")
  YX<-as.data.frame(YX)
  try(model01<-gls(Y~X,data=YX, correlation=bm))
  if(exists("model01")){
    summarymodel<-summary(model01)
    p.value.array_4th[snpIndex]<-summarymodel$tTable[,4][2]
  }else{p.value.array_4th[snpIndex]<-NA}
  rm(model01)
}
plot(-log10(p.value.array_4th))
save.image("drs4th.RData")

jpeg("drs4th.jpeg")
plot(-log10(p.value.array_4th))
dev.off()
