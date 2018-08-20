# http://www.gwaspi.org/?page_id=671
# To install core packages:
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# To install specific packages
# biocLite("snpStats")
#Read a PLINK binary data file as a SnpMatrix
#Description
#The package PLINK saves genome-wide association data in groups of three files, with the extensions .bed, .bim, and .fam. This function reads these files and creates an object of class "SnpMatrix"
rm(list=ls())
setwd("~/Documents/drosophila")
#install.packages("geiger")
library(snpStats)
library(nlme)
library(geiger)
library(ape)
library(phytools)
source("~/Dropbox/YiWeiHsu/Rcode/main.code/DV_test.r")
source("~/Dropbox/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
#snpmtx<-read.plink("dgrp2.bed","dgrp2.bim.txt","dgrp2.fam.txt", na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = NULL)
#snpgeno<-as.data.frame(snpmtx$genotypes)
#snpmap<-snpmtx$map
#rm(snpmtx)

load("snpdata.RData")
dim(snpgeno) #  4 million genome x 205  flies 

# chr1g<-snpgeno[,1:1013829]
# chr2g<-snpgeno[,1013830:1838775]
# chr3g<-snpgeno[,1838776:2840329]
# chr4g<-snpgeno[,2840330:3827021]
# chr5g<-snpgeno[,3827022:4427847]
chr6g<-snpgeno[,4427848:4438427]

dim(snpmap)  #  6 variables x 4 million genome   
# chr1m<-snpmap[snpmap$chromosome==1,]
# chr2m<-snpmap[snpmap$chromosome==2,]
# chr3m<-snpmap[snpmap$chromosome==3,]
# chr4m<-snpmap[snpmap$chromosome==4,]
# chr5m<-snpmap[snpmap$chromosome==5,]
chr6m<-snpmap[snpmap$chromosome==6,]

#we want to get the chromosome 6 for the flies 
#do chromosome 6

dim(chr6g)
dim(chr6m)

chr6gp<-matrix(as.numeric(t(chr6g)),ncol=205)
chr6gp<-t(chr6gp)
rownames(chr6gp)<-rownames(chr6g)
colnames(chr6gp)<-chr6m$position
table(chr6gp)
chr6gp[chr6gp==0]<-1
chr6gp[chr6gp==3]<-0
table(chr6gp)
dim(chr6gp)
#install.packages("ape")
#install.packages("phytools")

write.nexus.data(chr6gp,file="drs6th.nex")
system("paup -n get.tree.snp.nex >/dev/null")
system("grep 'PAUP_1' snp.tre | sed 's/^.*: //' > snp.1.tre ")
temp.tree<-read.table(file="snp.1.tre")
tree<-unlist(temp.tree[5])
tree<-as.character(tree)
write(tree,file="tree.txt")
tree<-read.newick(file="tree.txt")
writeNexus(tree,file="drs6th2.nex")
drs6th<-read.nexus("drs6th2.nex")

plot(tree)
trait<-read.csv("/Users/yiweihsu/Documents/drosophila/MackayStarvation/starvation.male.csv",header=FALSE)
head(trait)
rownames(trait)<-trait$V1
head(trait)
obj<-name.check(drs6th,trait)
obj
dim(trait)
drs6th.trait<-drop.tip(drs6th,obj$tree_not_data)
name.check(drs6th.trait,trait)
name.check(names(chr6gp[,1]),trait)
names(chr6gp[,1])

#chr6gpm<-chr6gp[match(drs6th.trait$tip.label, rownames(chr6gp) ),]

del1<-which(rownames(chr6gp)=="line_256")
del2<-which(rownames(chr6gp)=="line_509")
chr6gprm<-chr6gp[-c(del1,del2),]
bm<-corBrownian(1,drs6th.trait)
p.value.array6th<-array(0,c(dim(chr6gprm)[2]))
for(snpIndex in 1:dim(chr6gprm)[2]){
  print(snpIndex)
  YX<-cbind(trait$V2,chr6gprm[,snpIndex])
  colnames(YX)<-c("Y","X")
  YX<-as.data.frame(YX)
  try(model01<-gls(Y~X,data=YX, correlation=bm))
  if(exists("model01")){
  summarymodel<-summary(model01)
  p.value.array6th[snpIndex]<-summarymodel$tTable[,4][2]
  }else{p.value.array6th[snpIndex]<-NA}
  rm(model01)
  }
save.image("the6th.RData")

plot(-log10(p.value.array6th))
