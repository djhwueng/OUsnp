setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/OUsnp2_J2018")
rm(list=ls())
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV_test.r")
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")

library(ape)
library(mvtnorm)
library(invgamma)
library(geiger)
#put OU likelihood
VaF<-function(alpha,V=V){
  A<-exp(-2*alpha*(1-V))
  B<-(1-exp(-2*alpha*V))/(2*alpha)
  return(A*B)
  }

#Here is the likelihood profille
# OUnegloglike<-function(alpha,k=k,trait=trait,tree=tree){
#   DV.data<-DV(k=k,tree=tree)
#   V<-DV.data$V
#   D<-DV.data$D
#   Va<-VaF(alpha,V=V)
#   Va.inv<-pseudoinverse(Va)
#   mu.hat<-pseudoinverse(t(D)%*%Va.inv%*%D)%*%t(D)%*%Va.inv%*%trait
#   sigma.sq.hat<-t(trait-D%*%mu.hat)%*%Va.inv%*%(trait-D%*%mu.hat)/length(trait)
#   NegLogML <- (Ntip(tree)/2)*log(2*pi)+(1/2)*t(trait- D%*%mu.hat)%*%Va.inv%*%(trait- D%*%mu.hat) + (1/2)*log(abs(det(Va)))
#   return(NegLogML)
#   }

d.mtx<-function(n,snp=NULL){#formula (3) Thompson 2016
  Z.mtx <- array(0,c(n,2*n))
  #diplid index
  for(d.index in 1:n){
    Z.mtx[d.index,((2*d.index-1):(2*d.index))]<-0.5
   }
   return(Z.mtx)
  }
print(d.mtx(5))




########### Genetic component
diploid.size<-15
tree.size<-2*diploid.size
k<-5
#tree<-rcoal(tree.size)
seqlen<-1000
ouparams<-c(5,5,40,40,90,80,100)
names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
tree.trait<-sim.snp.tree.yg(taxa.size=tree.size,seqlen=seqlen,ouparams=ouparams)
tree<-tree.trait$tree
Zg<-matrix(tree.trait$trait$X,ncol=1) #here is Zg

load("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/BMsnp/test1.RData")
k<-5


alpha <- 0.5
DV.data<-DV(k=k,tree=tree)
V<-DV.data$V
Va<-VaF(alpha,V=V)
print(Va)
dim(Va)
D<-DV.data$D
mu<-array(1,c(k,1))
Z<-d.mtx(n=diploid.size)
V.yg<-Z%*%Va%*%t(Z)
print(V.yg)
E.yg<-Z%*%D%*%mu
sim.yg<-rmvnorm(n=1,mean=E.yg,sigma=V.yg)
sim.yg

####################### Environment component
nu.sq<-1
I<-diag(1,c(diploid.size))
r<-1
beta.vec<-array(c(1,2),c(r+1,1))
X.mtx<-array(0,c(diploid.size,r+1))
X.mtx[,1]<-1
X.mtx[,-1]<-rnorm(diploid.size*r)
print(X.mtx)
V.ye<-nu.sq*I
print(V.ye)
E.ye<-X.mtx%*%beta.vec
print(E.ye)
sim.ye<-rmvnorm(n=1,mean=E.ye,sigma=V.ye)

rho<-0.5
sim.y<-rho*sim.yg+(1-rho)*sim.ye
print(sim.y)


#shall add prior function posterior function proposal
# alpha.prior<-function(alpha){
#   return(dexp(alpha,rate=1,log=T))
#   }
# print(alpha.prior(1))
#as well as metropolis-hastings for alpha


likelihood<-function(model.params,Ye=Ye,Yg=Yg,r=r, k=k,rho=rho,tree=tree){
  alpha<-model.params[1]
  beta<-model.params[2:(r+2)]
  sigma.sq<-model.params[(r+3):(r+3)]
  nu.sq<-model.params[(r+4):(r+4)]
  mu<-model.params[(r+4+1):(r+4+k)]

  X<- matrix(runif(n=diploid.size,min=-5,max=5),ncol=1) #shall it be fixed all time ?
  one<-array(1,c(diploid.size,1))
  des.X<-cbind(one,X)
  #Ye.like<-dmvnorm( c(Ye) ,mean= des.X%*%beta, sigma= nu.sq*diag(1,c(length(Ye),length(Ye))),log=TRUE)

  Z<-d.mtx(n=length(Yg))
  DV.data<-DV(k=k,tree=tree)
  D<-DV.data$D
  V<-DV.data$V
  Va<-VaF(alpha,V=V)
  #Yg.like <- dmvnorm( c(Yg), mean=Z%*%D%*%mu, sigma=sigma.sq*pseudoinverse(Z%*%Va%*%t(Z)),log=TRUE)
  
  Y.loglike<-dmvnorm( c(Yg)+c(Ye), mean = des.X%*%beta + Z%*%D%*%mu ,  sigma = nu.sq*diag(1,c(length(Ye),length(Ye))) + sigma.sq*pseudoinverse(Z%*%Va%*%t(Z)) , log = TRUE)
  return( Y.loglike )
  }

#prior for gibbs
prior<-function(model.params,r=r,k=k,rho=rho){
  alpha<-model.params[1]
  beta<-model.params[2:(r+2)]
  sigma.sq<-model.params[(r+3):(r+3)]
  nu.sq<-model.params[(r+4):(r+4)]
  mu<-model.params[(r+4+1):(r+4+k)]

  alpha.prior<-dexp(alpha,rate=1,log=TRUE)

  u.sq<-100
  u.sq.I<-u.sq*diag(1,c(r+1,r+1))
  beta.prior<-dmvnorm(beta, mean=beta0, sigma=u.sq.I,log=TRUE) #PUT MVNNORM HERE

  #rho<-0.5 #will do a function input later
  a<-6
  b<-200*rho^2
  sigma.sq.prior <- dinvgamma(sigma.sq,shape=a,scale=b,log=TRUE)

  c<-6
  d<-200*(1-rho)^2
  nu.sq.prior<-dinvgamma(nu.sq,shape=c,scale=d,log=TRUE)

  mu0<-array(0,c(k,1))
  w.sq<-100
  w.sq.I<-w.sq*diag(1,c(k,k))
  mu.prior<-dmvnorm(mu,mean=mu0,sigma=w.sq.I,log=TRUE)

  return(alpha.prior+beta.prior+sigma.sq.prior+nu.sq.prior+mu.prior)
  }

posterior<-function(model.params,Yg=Yg,Ye=Ye,r=r,k=k,rho=rho,tree=tree){
    return(likelihood(model.params,Yg=Yg,Ye=Ye,r=r,k=k,rho=rho,tree=tree) + prior(model.params,r=r,k=k,rho=rho))
    }

alpha<-rexp(1)

beta0<- array(0,c(r+1,1)) #(14)
u.sq<- 100
u.sq.I<-u.sq*diag(1,c(r+1,r+1))
beta<-matrix(rmvnorm(n=1,mean=beta0,sigma=u.sq.I),ncol=1)

mu0<- array(0,c(k,1))  #(15)
w.sq<-100
w.sq.I<- w.sq*diag(1,c(k,k))
mu<-matrix(rmvnorm(n=1,mean=mu0,sigma=w.sq.I),ncol=1)

rho<-0.5

a<-6  #  (16)
b<-200*rho^2
sigma.sq <- rinvgamma(1,shape=a,scale=b)

c<-6 #c (17)
d<-200*(1-rho)^2 #d
nu.sq<-rinvgamma(1,shape=c,scale=d)

model.params<-c(alpha,beta, sigma.sq, nu.sq, mu)
names(model.params)<-c( "alpha",paste("b", 0:(length(beta)-1),sep="" ), "sig.sq","nu.sq",paste("mu",1:length(mu),sep=""))

#Ye<- rnorm(diploid.size)#later it shall change to (1-rho)*Ze
X<- matrix(runif(n=diploid.size,min=-5,max=5),ncol=1)
Ze<-matrix(rmvnorm(n=1,mean=20+4*X,sigma=40*diag(1,c(diploid.size, diploid.size))),ncol=1)
Ye<-(1-rho)*Ze

Zg<-Z%*%matrix(tree.trait$trait$X,ncol=1) #here is Zg
Zg
Yg<-rho*Zg
#Yg<- rnorm(diploid.size)#change to ms, rho*Zg  # we shall call Yg using ms.r

print(likelihood(model.params, Ye=Ye, Yg=Yg ,r=r, k=k,rho=rho, tree=tree))
print(prior(model.params,r=1,k=3, rho=0.5))
print(posterior(model.params, Ye=Ye, Yg=Yg ,r=r, k=k,rho=rho , tree=tree))

Y <- Yg+Ye

alpha.proposal<-function(alpha){
  # we can use MLE estimate and sd from geiger packages
  #?fitContinuous
  #trait<-matrix(trait,ncol=1)
  #tree<- rcoal(length(trait))#It would be better if we use the tree we construct. #compute.brlen(stree(length(trait),type="left"))
  #rownames(trait)<-tree$tip.label
  #ou.geiger<-fitContinuous(phy=tree,dat=trait, model="OU")
  #tree<-stree(length(trait),type="star")#can do this, need to fix DV function
  #ou.result<-optimize(OUnegloglike, lower=0,upper=10  ,k=k, trait=trait, tree=tree)
  #alpha
  #use star tree
  #return(rnorm(1,mean= ou.geiger$opt$alpha, sd=ou.geiger$opt$alpha/3 ) )
  return(rnorm(1, mean=alpha, sd = alpha/4))
  }
print(alpha.proposal(1))
#mh (alpha) within gibb
#Do

run_metropolis_MCMC_within_Gibbs<-function(startvalue=startvalue,iterations=interations,Yg=Yg,Ye=Ye,r=r,k=k, rho=rho,tree=tree, a=a,b=b,c=c,d=d){
  #iterations=50000
  chain=array(0,c(iterations+1,1+(r+1)+1+1+k)) #alpha, beta, sigma.sq, nu.sq, mu
  colnames(chain)<-c( "alpha",paste("beta",0:r,sep=""), "sigma.sq","nu.sq",paste("mu",1:k,sep=""))
  chain[1,]=startvalue
  diploid.size<-length(Yg)
  Z<-d.mtx(n=length(Yg))
  DV.data<-DV(k=k,tree=tree) #more care need on the DV function
  D<-DV.data$D
  V<-DV.data$V
  Va<- VaF(alpha,V=V)
  ZD<-Z%*%D
  X<-matrix(runif(n=diploid.size, min=5, max=5), ncol=1)
  one<-array(1,c(diploid.size,1))
  des.X<-cbind(one,X)
  tdes.X.des.X<-t(des.X)%*%des.X

  Ze<-matrix(rmvnorm(n=1,mean=20+4*X,sigma=40*diag(1,c(diploid.size, diploid.size))),ncol=1)
  Ye<-(1-rho)*Ze

  for(i in 2:iterations){

    if(i%%500==0)print(i)
    alpha<-chain[i-1,1] #need to check
    Va<-VaF(alpha,V=V)
    ZVatZ.inv<-pseudoinverse(Z%*%Va%*%t(Z))

    Y <- Yg+Ye

    #Gibb starts
    sigma.sq.shape <- a + diploid.size/2
    sigma.sq.scale <- b + 0.5*t(Yg-ZD%*%mu)%*%ZVatZ.inv%*%(Yg-ZD%*%mu)
    sigma.sq <- rinvgamma(1,shape=sigma.sq.shape, scale = sigma.sq.scale)

    nu.sq.shape<-c + diploid.size/2
    nu.sq.scale<-d + 0.5*t(Y-des.X%*%beta-Yg)%*%(Y-des.X%*%beta-Yg)
    nu.sq<-rinvgamma(1,shape=nu.sq.shape,scale=nu.sq.scale)

    V.yg<-pseudoinverse(ZVatZ.inv/sigma.sq  + diag(1,c(diploid.size,diploid.size))/nu.sq)
    V.yg.stuff<-V.yg%*%((Y-des.X%*%beta)/nu.sq  + ZVatZ.inv%*%ZD%*%mu/sigma.sq)
    Yg<-matrix(rmvnorm(n=1,mean=V.yg.stuff,sigma=V.yg),ncol=1)

    mu0<-array(0,c(k,1))
    V.mu <- pseudoinverse( t(ZD)%*%ZVatZ.inv%*%ZD/sigma.sq + diag(1,c(k,k))/w.sq )
    V.mu.stuff<-V.mu%*%(  t(ZD)%*%ZVatZ.inv%*%Yg /sigma.sq + mu0/w.sq  )
    mu<-matrix(rmvnorm(1,mean=V.mu.stuff,sigma=V.mu),ncol=1)

    u.sq<-100
    V.beta<-pseudoinverse(tdes.X.des.X/nu.sq + diag(1,dim(tdes.X.des.X))/u.sq )
    beta0<-array(0,c(r+1,1))
    V.beta.stuff<-V.beta%*%( t(des.X)%*%(Y-Yg) /nu.sq + beta0/u.sq )
    beta<-matrix(rmvnorm(n=1,mean=V.beta.stuff, sigma=V.beta), ncol=1)

    chain[i,2:length(startvalue)]<-c(beta, sigma.sq,nu.sq,mu)
    #Gibb ends

    #mh starts
    alpha.pps<-alpha.proposal(chain[i-1,1]) #slower
    proposal<-c(alpha.pps, chain[i,2:length(startvalue)])

#    likelihood(proposal,Yg=Yg,Ye=Ye,r=r,k=k,rho=rho,tree=tree)
#    model.params<-proposal


    probab=exp(posterior(proposal,Yg=Yg,Ye=Ye,r=r,k=k, rho=rho,tree=tree) - posterior(chain[i-1,],Yg=Yg,Ye=Ye,r=r,k=k, rho=rho,tree=tree))
      if(runif(1)<probab){
        chain[i,]=c(alpha.pps, chain[i,2:length(startvalue)])
        }else{
          chain[i,]=c(chain[i-1,1], chain[i,2:length(startvalue)] )
        }
    #mh ends
      }#end for loop
  return(chain)
  }

#alpha            b0            b1        sig.sq         nu.sq           mu1           mu2       mu3
# 0.391505446   0.394397955  -8.062240686   0.004365570   0.002004884 -12.779226547  -5.700671144  -10.547468962

startvalue= c(rexp(1), rnorm(r+1), rexp(1), rexp(1), rnorm(k) )

chain = run_metropolis_MCMC_within_Gibbs(startvalue=startvalue,iterations=50000,Yg=Yg,Ye=Ye,r=r,k=k, rho=rho,tree=tree, a=a,b=b,c=c,d=d)
burnIn =10000
acceptance = 1- mean(duplicated(chain[-(1:burnIn),] ))
print(acceptance)

par(mfrow=c(2,4))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of alpha", xlab= "True value = red line"  )
abline(v= mean(chain[-(1:burnIn),1]))
abline(v=alpha,col="red")

hist(chain[-(1:burnIn),1+1+r],nclass=30, main="Posterior of beta1", xlab= "True value = red line"  )
abline(v= mean(chain[-(1:burnIn),1+1+r]))
abline(v=beta[2],col="red")

hist(chain[-(1:burnIn),1+1+r+1],nclass=30, main="Posterior of sigma.sq", xlab= "True value = red line"  )
abline(v= mean(chain[-(1:burnIn),1+1+r+1]))
abline(v=sigma.sq,col="red")

hist(chain[-(1:burnIn),1+1+r+1+1],nclass=30, main="Posterior of nu.sq", xlab= "True value = red line"  )
abline(v= mean(chain[-(1:burnIn),1+1+r+1+1]))
abline(v=nu.sq,col="red")


plot(chain[-(1:burnIn),1],type="l",xlab="True value of red line", main="alpha")
abline(h=alpha, col="red")

plot(chain[-(1:burnIn),1+1+r],type="l",xlab="True value of red line", main="beta1")
abline(h=beta[2], col="red")

plot(chain[-(1:burnIn),1+1+r+1],type="l",xlab="True value of red line", main="sigma.sq")
abline(h=sigma.sq, col="red")

plot(chain[-(1:burnIn),1+1+r+1+1],type="l",xlab="True value of red line", main="nu.sq")
abline(h=nu.sq, col="red")
system("pwd")
save.image("test1.RData")
