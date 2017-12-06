rm(list=ls())
library(phytools)
library(phyclust)
library(phangorn)

ntip<-5
phy<-pbtree(n=ntip)
anc<-phy$edge[,1]
des<-phy$edge[,2]
N<-dim(phy$edge)[1]
bm.sim.trait <- rep(0, ntip + phy$Nnode)
ou.sim.trait <- rep(0, ntip + phy$Nnode)
ROOT<-ntip+1
bm.sim.trait[ROOT]<- 0 #start from the mean, maybe average trait value
ou.sim.trait[ROOT]<- 0 #it could be the formula mu+(root-mu)exp(-2alpha t)
sigma.sq<-1
alpha<-0.00001

for(i in N:1){
#  i<-N-9
#  print(des[i])
#  print(anc[i])
  bm.sim.trait[des[i]]<- bm.sim.trait[anc[i]]+ sigma.sq*rnorm(n=1,mean=0,sd=sqrt(phy$edge.length[i]))
  ou.mean<-ou.sim.trait[ROOT] + (ou.sim.trait[anc[i]]-ou.sim.trait[ROOT])*exp(-alpha*phy$edge.length[i]) 
  ou.var<- sigma.sq/(2*alpha)*(1-exp(-2*alpha*phy$edge.length[i]))
  #print(c(phy$edge.length[i],ou.var))
  ou.sim.trait[des[i]]<-ou.sim.trait[anc[i]]+ rnorm(n=1,mean=ou.mean,sd=sqrt(ou.var))
  }

print(bm.sim.trait)
print(ou.sim.trait)
#   
#   #    x
#       #print( c(anc[i],des[i] ))
#   #print(  x[des[i]])
#   }
# 

#sim dna
#https://snoweye.github.io/phyclust/document/html/ms.html
# ret.ms<-ms(nsam=3,nreps=2,opts="-T -G 0.1")
# tree.anc<-read.tree(text=ret.ms[3])
# tree.anc
# tree.anc$tip.label<-paste("a",1:Ntip(tree.anc),sep="")
# plot(tree.anc)
# 
# dna<-seqgen(opts="-mHKY -l10 -s0.2",newick.tree=ret.ms[3])
# print(dna)

#sim snp
?gen.seq.SNP()
ret.ms<-ms(nsam=5,nreps=1,opts="-T")
ret.ms
tree.ms<-read.tree(text=ret.ms[3])
tree.ms
anc.SNP<-rep(0:1,6)
anc.SNP
pi.SNP<-c(0.4,0.6)
pi.SNP
L<-length(anc.SNP)
L
set.seed(1234)
paste(sid2snp(anc.SNP),collapse = "")
SNP.1<-gen.seq.SNP(tree.ms,pi.SNP,L,anc.seq=anc.SNP)
SNP.1
SNP.2<-gen.seq.SNP(tree.ms,pi.SNP,L,rate.scale=5,anc.seq=anc.SNP)
SNP.2
