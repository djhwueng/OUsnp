setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/BMsnp")
rm(list=ls())
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV.r")
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
VaF<-function(alpha,V=V){
  A<-exp(-2*alpha*(1-V))
  B<- (1-exp(-2*alpha*V))/(2*alpha)
  return(A*B)
  }

#for OU we need optimization here
OUnegloglike<-function(alpha, k=k, trait=trait, tree=tree){
    DV.data<-DV(k=k,tree=tree)
    V<-DV.data$V
    D<-DV.data$D
    Va<-VaF(alpha,V=V)
    Va.inv<-pseudoinverse(Va)
    mu.hat<-pseudoinverse(t(D)%*%Va.inv%*%D)%*%t(D)%*%Va.inv%*%trait
    sigma.sq.hat<- t(trait - D%*%mu.hat)%*%Va.inv%*%(trait - D%*%mu.hat)/ntip
#    print(mu.hat)
#    print(sigma.sq.hat)
    NegLogML <- (Ntip(tree)/2)*log(2*pi)+(1/2)*t(trait- D%*%mu.hat)%*%Va.inv%*%(trait- D%*%mu.hat) + (1/2)*log(abs(det(Va)))
    return(NegLogML)
    }

MLEsOU <- function(alpha, k=k, trait=trait, tree=tree){
  DV.data<-DV(k=k,tree=tree)
  V<-DV.data$V
  D<-DV.data$D
  Va<-VaF(alpha,V=V)
  Va.inv<-pseudoinverse(Va)
  mu.hat<-pseudoinverse(t(D)%*%Va.inv%*%D)%*%t(D)%*%Va.inv%*%trait
  sigma.sq.hat<- t(trait - D%*%mu.hat)%*%Va.inv%*%(trait - D%*%mu.hat)/ntip
  return(list(mu.hat=mu.hat,sigma.sq.hat=sigma.sq.hat,alpha.hat=alpha))
  }

LSS.ou<-function(k,trait=trait,tree=tree){
  n<-length(trait)
  ou.result<-optimize(OUnegloglike, lower=0,upper=10  ,k=k, trait=trait, tree=tree)
  #ou.mle<-MLEsOU(c(ou.result$objective), k=k, trait=trait, tree=tree)
  loglik<- -ou.result$minimum
  lss<- 2*loglik - k*n
  return(lss)
  }

permute.LSS.ou<- function(k,trait=trait,tree=tree,sims=sims){
  lss.array<-array(0,c(sims))
  for(sim.index in 1:sims){
    print(sim.index)
    permute.trait<-sample(trait)
    permute.lss.ou <- LSS.ou(k,trait=permute.trait,tree=tree)
    lss.array[sim.index]<-permute.lss.ou
    }
  return(lss.array)
  }


#k<-3
k.array<-1:5
ntip<-10
#tree<-pbtree(n=ntip)
#trait<-rnorm(ntip)

taxa.size<-10
seqlen<-1000
ouparams<-c(5,5,40,40,90,80,100)
names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
tree.trait<-sim.snp.tree.yg(taxa.size=taxa.size,seqlen=seqlen,ouparams=ouparams)
tree<-tree.trait$tree
trait<-matrix(tree.trait$trait,ncol=1)


obs.LSS.ou.array<- array(0,c(length(k.array)))
for(k in 1:length(k.array)){
  obs.LSS.ou.array[k] <- LSS.ou(k=k, trait=trait, tree=tree)
  }
print(obs.LSS.ou.array)
max.LSS.ou<-max(obs.LSS.ou.array)
print(max.LSS.ou)
#alpha<- c(obs.LSS.ou.array$objective)

sims<-100
sim.LSS.ou.array<-array(0,c(k,sims))
for(simIndex in 1:100){
  trait<-sample(trait)
  for(k in 1:length(k.array)){
    sim.LSS.ou.array[k,simIndex]<-LSS.ou(k=k, trait=trait, tree=tree)
    }
 }
sim.max.LSS.ou<-apply(sim.LSS.ou.array, 2 ,max)
sim.max.LSS.ou
#ou.result<-optimize(OUnegloglike, lower=0,upper=10  ,k=k, trait=trait, tree=tree)
#print(ou.result)
# alpha<- c(ou.result$objective)
#print(alpha)
#print(OUnegloglike(alpha, k=k, trait=trait, tree=tree))
# a bit wierd here for the likelihood, expect to check later.
# print(MLEsOU(alpha, k=k, trait=trait, tree=tree))
# obs.LSS.ou<-LSS.ou(k,trait = trait , tree = tree)
# print(LSS.ou(k,trait = trait , tree = tree))

# sims<-100
# permute.LSS.ou.array<-permute.LSS.ou(k,trait=trait,tree=tree,sims=sims)
#the p.value for detection at each locus is the proportion of
#data sets scoring higher than the observed data set at each particular locus.
p.value<-sum(sim.LSS.ou.array>c(max.LSS.ou))
print(p.value)
