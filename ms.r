rm(list=ls())
library(phytools)
library(gap)
library(phyclust)
library(OUwie)
library(ape)
# ?BMhyb
system("pwd")
# setwd("~/Dropbox/YiWeiHsu/ms")
setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/ms")

taxa.size<-5
sim.snp.data<-paste("ms ",taxa.size, " 1 -T -s 100 > ms.out" ,sep = "")

#system("ms 10 1 -T -s 100  > ms.out")
system(sim.snp.data)
system("less ms.out")
msout1<- read.ms.output("ms.out")
# msout1
ls(msout1)
snp<-matrix(unlist(msout1$gametes),nrow=taxa.size)
snp
write.nexus.data(snp,file="snp100.nex",format="dna")
system("less snp100.nex")
system("sed 's/DNA/standard/g' snp100.nex > snp100_1.nex")
system("less snp100_1.nex")
system("less test_gettree.snp.2.nex")
system("paup -n test_gettree.snp.2.nex > /dev/null")
system("grep 'PAUP_1' snp100_2.tre | sed 's/^.*: //' >2.tre")
system("less 2.tre")
temp.tree<-read.table(file="2.tre")
temp.tree
typeof(temp.tree[5])
tree<-unlist(temp.tree[5])
tree<-as.character(tree)
#tree<-"(taxon_1:16,(((taxon_2:26,taxon_4:30):14,taxon_5:7):17,(((taxon_3:17,taxon_7:36):15,taxon_6:8):27,taxon_8:15):15):10);"
write(tree, file="tree.txt")
system("less tree.txt")
tree<-read.newick(file="tree.txt")
tree
plot(tree)
vcv(tree)
# is.ultrametric(tree)
tree$root.edge<-0

alpha=c(5,5)
sigma.sq=c(10,20)
theta0=0.0
theta=c(80,100)
tree.for.simmap<-sim.history(tree,Q=matrix(c(-1,1,1,-1),2,2))#Juke-Cantor
plot(tree.for.simmap)
sim.data<-OUwie.sim(tree.for.simmap,simmap.tree=TRUE,scaleHeight=FALSE,root.age=0,
                  alpha=alpha,sigma.sq=sigma.sq,theta0=theta0,theta=theta)
print(sim.data)


