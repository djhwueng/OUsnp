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

chr5g<-snpgeno[,3827022:4427847]
chr5m<-snpmap[snpmap$chromosome==5,]
chr5gp<-matrix(as.numeric(t(chr5g)),ncol=205)
chr5gp<-t(chr5gp)
rownames(chr5gp)<-rownames(chr5g)
colnames(chr5gp)<-chr5m$position
table(chr5gp)
chr5gp[chr5gp==0]<-1
chr5gp[chr5gp==3]<-0

write.nexus.data(chr5gp,file="drs5.nex")
system("sed 's/DNA/standard/g' drs5.nex >drs5.1.nex")
system("paup -n get.tree.snp.drs5.nex >/dev/null")
system("grep 'PAUP_1' drs5.tre | sed 's/^.*: //' > drs5.1.tre ")
temp.tree<-read.table(file="drs5.1.tre")
tree<-unlist(temp.tree[5])
tree<-as.character(tree)
write(tree,file="treedrs5.txt")
tree<-read.newick(file="treedrs5.txt")
writeNexus(tree,file="drs5.nex")


drs5<-read.nexus("drs5.nex")
plot(drs5)
#trait<-read.csv("/Users/yiweihsu/Documents/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)
trait<-read.csv("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)
head(trait)
rownames(trait)<-trait$V1
head(trait)
obj<-name.check(drs5,trait)
obj

dim(trait)
drs5.trait<-drop.tip(drs5,obj$tree_not_data)
name.check(drs5.trait,trait)


del1<-which(rownames(chr5gp)=="line_256")
del2<-which(rownames(chr5gp)=="line_509")
delset<-c(del1,del2)

chr5gprm<-chr5gp[-delset,]
bm<-corBrownian(1,drs5.trait)
p.value.array_5th<-array(0,c(dim(chr5gprm)[2]))
for(snpIndex in 1:dim(chr5gprm)[2]){
  print(snpIndex)
  YX<-cbind(trait$V2,chr5gprm[,snpIndex])
  colnames(YX)<-c("Y","X")
  YX<-as.data.frame(YX)
  try(model01<-gls(Y~X,data=YX, correlation=bm))
  if(exists("model01")){
    summarymodel<-summary(model01)
    p.value.array_5th[snpIndex]<-summarymodel$tTable[,4][2]
  }else{p.value.array_5th[snpIndex]<-NA}
  rm(model01)
}
plot(-log10(p.value.array_5th))
save.image("drs5th.RData")

jpeg("drs5th.jpeg")
plot(-log10(p.value.array_5th))
dev.off()
