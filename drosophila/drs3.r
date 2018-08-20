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
setwd("~/Documents/drosophila")
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
source("~/Dropbox/YiWeiHsu/Rcode/main.code/DV_test.r")
source("~/Dropbox/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
# snpmtx<-read.plink("dgrp2.bed","dgrp2.bim.txt","dgrp2.fam.txt", na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = NULL)
# snpgeno<-as.data.frame(snpmtx$genotypes)
# snpmap<-snpmtx$map
# rm(snpmtx)

load("snpdata.RData")
dim(snpgeno) #  4 million genome x 205  flies 

chr3g<-snpgeno[,1838776:2840329]
chr3m<-snpmap[snpmap$chromosome==3,]
chr3gp<-matrix(as.numeric(t(chr3g)),ncol=205)
chr3gp<-t(chr3gp)
rownames(chr3gp)<-rownames(chr3g)
colnames(chr3gp)<-chr3m$position
table(chr3gp)
chr3gp[chr3gp==0]<-1
chr3gp[chr3gp==3]<-0

write.nexus.data(chr3gp,file="drs3.nex")
system("sed 's/DNA/standard/g' drs3.nex >drs3.1.nex")
system("paup -n get.tree.snp.drs3.nex >/dev/null")
system("grep 'PAUP_1' drs3.tre | sed 's/^.*: //' > drs3.1.tre ")
temp.tree<-read.table(file="drs3.1.tre")
tree<-unlist(temp.tree[5])
tree<-as.character(tree)
write(tree,file="treedrs3.txt")
tree<-read.newick(file="treedrs3.txt")
writeNexus(tree,file="drs3.nex")
drs3<-read.nexus("drs3.nex")
plot(drs3)
trait<-read.csv("/Users/yiweihsu/Documents/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)
head(trait)
rownames(trait)<-trait$V1
head(trait)
obj<-name.check(drs3,trait)
obj
dim(trait)
drs3.trait<-drop.tip(drs3,obj$tree_not_data)
name.check(drs3.trait,trait)


del1<-which(rownames(chr3gp)=="line_256")
del2<-which(rownames(chr3gp)=="line_509")
delset<-c(del1,del2)

chr3gprm<-chr3gp[-delset,]
bm<-corBrownian(1,drs3.trait)
p.value.array_3th<-array(0,c(dim(chr3gprm)[2]))
for(snpIndex in 1:dim(chr3gprm)[2]){
  print(snpIndex)
  YX<-cbind(trait$V2,chr3gprm[,snpIndex])
  colnames(YX)<-c("Y","X")
  YX<-as.data.frame(YX)
  try(model01<-gls(Y~X,data=YX, correlation=bm))
  if(exists("model01")){
    summarymodel<-summary(model01)
    p.value.array_3th[snpIndex]<-summarymodel$tTable[,4][2]
  }else{p.value.array_3th[snpIndex]<-NA}
  rm(model01)
}
plot(-log10(p.value.array_3th))
save.image("drs3th.RData")

jpeg("drs3th.jpeg")
plot(-log10(p.value.array_3th))
dev.off()
