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

chr2g<-snpgeno[,1013830:1838775]
chr2m<-snpmap[snpmap$chromosome==2,]
chr2gp<-matrix(as.numeric(t(chr2g)),ncol=205)
chr2gp<-t(chr2gp)
rownames(chr2gp)<-rownames(chr2g)
colnames(chr2gp)<-chr2m$position
table(chr2gp)
chr2gp[chr2gp==0]<-1
chr2gp[chr2gp==3]<-0

write.nexus.data(chr2gp,file="drs2.nex")
system("sed 's/DNA/standard/g' drs2.nex >drs2.1.nex")
system("paup -n get.tree.snp.drs2.nex >/dev/null")
system("grep 'PAUP_1' drs2.tre | sed 's/^.*: //' > drs2.1.tre ")
temp.tree<-read.table(file="drs2.1.tre")
tree<-unlist(temp.tree[5])
tree<-as.character(tree)
write(tree,file="treedrs2.txt")
tree<-read.newick(file="treedrs2.txt")
writeNexus(tree,file="drs2.nex")



drs2<-read.nexus("drs2.nex")
plot(drs2)
trait<-read.csv("/Users/yiweihsu/Documents/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)
head(trait)
rownames(trait)<-trait$V1
head(trait)
obj<-name.check(drs2,trait)
obj
dim(trait)
drs2.trait<-drop.tip(drs2,obj$tree_not_data)
name.check(drs2.trait,trait)


del1<-which(rownames(chr2gp)=="line_256")
del2<-which(rownames(chr2gp)=="line_509")
delset<-c(del1,del2)

chr2gprm<-chr2gp[-delset,]
bm<-corBrownian(1,drs2.trait)
p.value.array_2th<-array(0,c(dim(chr2gprm)[2]))
for(snpIndex in 1:dim(chr2gprm)[2]){
  print(snpIndex)
  YX<-cbind(trait$V2,chr2gprm[,snpIndex])
  colnames(YX)<-c("Y","X")
  YX<-as.data.frame(YX)
  try(model01<-gls(Y~X,data=YX, correlation=bm))
  if(exists("model01")){
    summarymodel<-summary(model01)
    p.value.array_2th[snpIndex]<-summarymodel$tTable[,4][2]
  }else{p.value.array_2th[snpIndex]<-NA}
  rm(model01)
}
plot(-log10(p.value.array_2th))
save.image("drs2th.RData")

jpeg("drs2th.jpeg")
plot(-log10(p.value.array_2th))
dev.off()
