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

chr1g<-snpgeno[,1:1013829]
chr1m<-snpmap[snpmap$chromosome==1,]
chr1gp<-matrix(as.numeric(t(chr1g)),ncol=205)
chr1gp<-t(chr1gp)
rownames(chr1gp)<-rownames(chr1g)
colnames(chr1gp)<-chr1m$position
table(chr1gp)
chr1gp[chr1gp==0]<-1
chr1gp[chr1gp==3]<-0

write.nexus.data(chr1gp,file="drs1.nex")
system("sed 's/DNA/standard/g' drs1.nex >drs1.1.nex")
system("paup -n get.tree.snp.drs1.nex >/dev/null")
system("grep 'PAUP_1' drs1.tre | sed 's/^.*: //' > drs1.1.tre ")
temp.tree<-read.table(file="drs1.1.tre")
tree<-unlist(temp.tree[5])
tree<-as.character(tree)
write(tree,file="treedrs1.txt")
tree<-read.newick(file="treedrs1.txt")
writeNexus(tree,file="drs1.nex")


drs1<-read.nexus("drs1.nex")
plot(drs1)
trait<-read.csv("/Users/yiweihsu/Documents/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)
head(trait)
rownames(trait)<-trait$V1
head(trait)
obj<-name.check(drs1,trait)
obj
dim(trait)
drs1.trait<-drop.tip(drs1,obj$tree_not_data)
name.check(drs1.trait,trait)


del1<-which(rownames(chr1gp)=="line_256")
del2<-which(rownames(chr1gp)=="line_509")
delset<-c(del1,del2)

chr1gprm<-chr1gp[-delset,]
bm<-corBrownian(1,drs1.trait)
p.value.array_1th<-array(0,c(dim(chr1gprm)[2]))
for(snpIndex in 1:dim(chr1gprm)[2]){
  print(snpIndex)
  YX<-cbind(trait$V2,chr1gprm[,snpIndex])
  colnames(YX)<-c("Y","X")
  YX<-as.data.frame(YX)
  try(model01<-gls(Y~X,data=YX, correlation=bm))
  if(exists("model01")){
    summarymodel<-summary(model01)
    p.value.array_1th[snpIndex]<-summarymodel$tTable[,4][2]
  }else{p.value.array_1th[snpIndex]<-NA}
  rm(model01)
}
plot(-log10(p.value.array_1th))
save.image("drs1th.RData")

jpeg("drs1th.jpeg")
plot(-log10(p.value.array_1th))
dev.off()
