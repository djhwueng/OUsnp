rm(list=ls())
setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/BMsnp")
# setwd("~/Dropbox/YiWeiHsu/Rcode/main.code/BMsnp")
 source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV_test.r")
# source("~/Dropbox/YiWeiHsu/Rcode/main.code/DV_test.r")
 source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
#source("~/Dropbox/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")


LSS.bm<-function(k=k,trait=trait,V=V,D=D){
  n<-length(trait)
  V.inv<-pseudoinverse(V)
  mu.hat<-pseudoinverse(t(D)%*%V.inv%*%D)%*%t(D)%*%V.inv%*%trait
  sigma.sq.hat<- t(trait - D%*%mu.hat)%*%V.inv%*%(trait - D%*%mu.hat)/n
  NegLogML <- (n/2)*log(2*pi)+(1/2)*t(trait- D%*%mu.hat)%*%V.inv%*%(trait- D%*%mu.hat) + (1/2)*log(abs(det(V))) #natural log multivariate noraml PDF取log
  loglik<- -NegLogML
  lss<- 2*loglik-k*log(n)
  return(lss)
  }

permute.LSS.bm<- function(k,trait=trait,tree=tree,sims=sims){
  lss.array<-array(0,c(sims))
  for(sim.index in 1:sims){
    #    print(sim.index)
    permute.trait<-sample(trait)
    permute.lss.bm <- LSS.bm(k,trait=permute.trait,tree=tree)
    lss.array[sim.index]<-permute.lss.bm
     }
  return(lss.array)
  }

#tree<-pbtree(n=5)#this is a sumulated tree
#plot(tree)
#nodelabels()
#trait<-rnorm(ntip) # this is a simulated trait
taxa.size<-30
seqlen<-1000
ouparams<-c(5,5,40,40,90,80,100)
names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
tree.trait<-sim.snp.tree.yg(taxa.size=taxa.size,seqlen=seqlen,ouparams=ouparams)
tree<-tree.trait$tree
plot(tree)
#par(mfrow=c(1,2))
#plot(tree)
#tree<-reorder(tree,"postorder")
#nodelabels()
tree<-chronoMPL(tree)#we
plot(tree)
#nodelabels()

tree$tip.label<-paste("t",tree$tip.label,sep="")
tree$tip.label


vcv.tree<-vcv(tree)
diag(vcv.tree)<-max(vcv.tree)
vcv.tree<-vcv.tree/max(vcv.tree)
tree<-vcv2phylo(vcv.tree)
plot(tree)
trait<-matrix(tree.trait$trait$X,ncol=1)
rownames(trait)<-tree.trait$trait$Genus_species
trait

#we would like to get D and V first and store it.

k.array <- 3:6
# for(kIndex in 1:length(k.array)){
#     k<-k.array[kIndex]
#     k<-3
#     #k<-4
#     #DV.data<-
#     print(k)
#     DV.data<-DV(k=k,tree=tree)
#     V<-DV.data$V
#     D<-DV.data$D
#     dim(V)
#     dim(D)
#     D
#     #print(dim(D)[2])
#     }


obs.LSS.bm.array <- array(0,c(length(k.array)))
for(kIndex in 1:length(k.array)){
  print(kIndex)
  DV.data<-DV(k=k.array[kIndex],tree=tree)
  V<-DV.data$V
  D<-DV.data$D
  print(dim(V))
  dim(D)
  obs.LSS.bm.array[kIndex] <- LSS.bm(k=k.array[kIndex], trait=trait, V=V,D=D)
  }

print(obs.LSS.bm.array)
plot(obs.LSS.bm.array)
max.LSS.bm<-max(obs.LSS.bm.array)
print(max.LSS.bm)
# DV.data<-DV(k=k,tree=tree)
# V<-DV.data$V
# D<-DV.data$D
# V.inv<-pseudoinverse(V)
# mu.hat<-pseudoinverse(t(D)%*%V.inv%*%D)%*%t(D)%*%V.inv%*%trait
# sigma.sq.hat<- t(trait - D%*%mu.hat)%*%V.inv%*%(trait - D%*%mu.hat)/ntip
# print(mu.hat)
# print(sigma.sq.hat)
#
#
# obs.LSS.bm<-LSS.bm(k,trait=trait,tree=tree)
# print(LSS.bm(k,trait=trait,tree=tree))
#
#

rm(list=ls())
load("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/BMsnp/test1.RData")

sims<-1000
sim.LSS.bm.array<-array(0,c(length(k.array),sims))
max.sim.LSS.bm.array<-array(0,c(sims))
for(simIndex in 1:sims){
   trait<-sample(trait)
   for(kIndex in 1:length(k.array)){
   sim.LSS.bm.array[kIndex,simIndex]<-LSS.bm(k=k.array[kIndex],trait=trait,V=V,D=D)
      }
   max.sim.LSS.bm.array[simIndex]<-max(sim.LSS.bm.array[,simIndex])
     }

 #permute.LSS.bm.array<-permute.LSS.bm(k,trait=trait,tree=tree,sims=sims)
 #the p.value for detection at each locus is the proportion of
 # #data sets scoring higher than the observed data set at each particular locus.
p.value<-sum(max.sim.LSS.bm.array>c(max.LSS.bm))
print(p.value)
#
# sample(1:5)

system("pwd")
save.image(file="test1.RData")
# 將原始資料丟進LSS裡面跑，出來的最大值就是我們分的群數，然後我們將分數重新做排列(sample)，重新排列完
# 在丟回LSS裡面做計算，看看兩次的分數，有差異的地方就是有顯著性的SNP
