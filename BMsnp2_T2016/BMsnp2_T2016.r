setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/BMsnp2_T2016")
rm(list=ls())
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV.r")
source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
library(ape)
library(mvtnorm)
library(invgamma)
#For sth diploid individual, the genetic component of the
#quantitative trait, Ygs = Ys, can be represented as the mean of two allelic
#trait values, Ys1 and Ys2. The covariance among any two allelic trait values
#among any two allelic trait values can be determined
#using the phylogenetic tree at the considered SNP.
d.mtx<-function(n,snp=NULL){ #formula (3) Thompson 2016
  Z.mtx<-array(0,c(n,2*n))
  #diploid index
  for(d.index in 1:n){
    Z.mtx[d.index,((2*d.index-1):(2*d.index))]<-0.5
    Z.mtx
    }
  return(Z.mtx)
  }

print(d.mtx(5))

diploid.size<-5 #laura did size of 50
#trait<-array(rnorm(diploid.size),c(diploid.size,1))

############## Genetic component
tree.size<-2*diploid.size
k<-3 #cluster, it is a fixed number in Thompson 2016. laura did 5
sigma.sq<-1
tree<-rcoal(tree.size)

seqlen<-1000
ouparams<-c(5,5,40,40,90,80,100)
names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
tree.trait<-sim.snp.tree.yg(taxa.size=tree.size,seqlen=seqlen,ouparams=ouparams)
tree<-tree.trait$tree
Zg<-matrix(tree.trait$trait,ncol=1) #here is Zg


DV.data<-DV(k=k,tree=tree)
V<-DV.data$V
D<-DV.data$D
mu<-array(1,c(k,1))
#V.inv<-pseudoinverse(V)
#mu.hat<-pseudoinverse(t(D)%*%V.inv%*%D)%*%t(D)%*%V.inv%*%trait
Z<-d.mtx(n=diploid.size)
V.yg<-sigma.sq*Z%*%V%*%t(Z)
print(V.yg)
E.yg<-Z%*%D%*%mu
print(E.yg)
#simulate yg
#https://cran.r-project.org/web/packages/mvtnorm/mvtnorm.pdf
sim.yg <- rmvnorm(n=1,mean=E.yg,sigma=V.yg)

############## Environment component
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
sim.ye <- rmvnorm(n=1,mean=E.ye,sigma=V.ye)

############ Trait
rho<-0.5
sim.y<- rho*sim.yg + (1-rho)*sim.ye
print(sim.y)

####ABOVE ARE THE CODE FOR CHECKING CODE

############## prior
beta0<- array(0,c(r+1,1)) #(14)
u.sq<- 100
u.sq.I<-u.sq*diag(1,c(r+1,r+1))
beta<-matrix(rmvnorm(n=1,mean=beta0,sigma=u.sq.I),ncol=1)

mu0<- array(0,c(k,1))  #(15)
w.sq<-100
w.sq.I<- w.sq*diag(1,c(k,k))
mu<-matrix(rmvnorm(n=1,mean=mu0,sigma=w.sq.I),ncol=1)

a<-6  #  (16)
b<-200*rho^2
sigma.sq <- rinvgamma(1,shape=a,scale=b)

c<-6 #c (17)
d<-200*(1-rho)^2 #d
nu.sq<-rinvgamma(1,shape=c,scale=d)



model.params<-c(beta, sigma.sq, nu.sq, mu)
names(model.params)<-c( paste("b", 0:(length(beta)-1),sep="" ), "sig.sq","nu.sq",paste("mu",1:length(mu),sep=""))


sim.chains<-50000
gibb.sample<-array(0,c(sim.chains,length(model.params)))
colnames(gibb.sample)<-c( paste("b", 0:(length(beta)-1),sep="" ), "sig.sq","nu.sq",paste("mu",1:length(mu),sep=""))
gibb.sample[1,]<-model.params

##############
X<- matrix(runif(n=diploid.size,min=-5,max=5),ncol=1)
Ze<-matrix(rmvnorm(n=1,mean=20+4*X,sigma=40*diag(1,c(diploid.size, diploid.size))),ncol=1)
Ye<-(1-rho)*Ze
#########################
#########################
Zg <- rnorm(diploid.size) #need use ms, tree, OUwie ... later
#########################
#########################
Yg<-rho*Zg

Y <- Yg+Ye
#starting point here from prior
#loop and array
ZVtZ.inv<-pseudoinverse(Z%*%V%*%t(Z))
one<-array(1,c(diploid.size,1))
des.X<-cbind(one,X)
tdes.X.des.X<- t(des.X)%*%des.X
ZD<-Z%*%D

for(chainIndex in 2:sim.chains){
    Y<-Yg+Ye
    if(chainIndex %%100 ==0){print(chainIndex)}
    sigma.sq.shape <- a + diploid.size/2
    sigma.sq.scale <- b + 0.5*t(Yg-ZD%*%mu)%*%ZVtZ.inv%*%(Yg-ZD%*%mu)
    sigma.sq<- rinvgamma(1, shape=sigma.sq.shape,scale=sigma.sq.scale)#(19)  #sig.sq | y,yg,beta,mu,nu.sq

    nu.sq.shape <- c + diploid.size/2
    nu.sq.scale <- d + 0.5*t(Y- des.X%*%beta-Yg)%*%(Y- des.X%*%beta-Yg)
    nu.sq<-rinvgamma(1,shape=nu.sq.shape,scale=nu.sq.scale) #(20)

    V.yg<-pseudoinverse(ZVtZ.inv/sigma.sq + diag(1,c( diploid.size,diploid.size))/nu.sq)
    V.yg.stuff<-V.yg%*%( (Y-des.X%*%beta)/nu.sq  +  ZVtZ.inv%*%ZD%*%mu/sigma.sq)
    Yg<-matrix(rmvnorm(n=1,mean=V.yg.stuff,sigma=V.yg),ncol=1) #(22)

    V.mu <- pseudoinverse(t(ZD)%*%ZVtZ.inv%*%ZD/sigma.sq + diag(1,c(k,k))/w.sq )
    V.mu.stuff <-  V.mu%*% (  t(ZD)%*%ZVtZ.inv%*%Yg  /sigma.sq   +    mu0/w.sq)
    mu<- matrix(rmvnorm(n=1,mean= V.mu.stuff, sigma=V.mu),ncol=1)#(23)
              #(24)
    V.beta<-pseudoinverse(  tdes.X.des.X /nu.sq + diag(1,dim(tdes.X.des.X) )/u.sq)
    V.beta.stuff<-V.beta%*%( t(des.X)%*%(Y-Yg)   /nu.sq +  beta0 /u.sq )
    beta <- matrix(rmvnorm(n=1, mean = V.beta.stuff, sigma=V.beta),ncol=1)    #(25)

    gibb.sample[chainIndex,]<-c(beta, sigma.sq, nu.sq, c(mu))
    }

burnIn=10000
#acceptance = 1- mean(duplicated(gibb.sample[-(1:burnIn),]))
#print(acceptance)
#rinvgamma(1,7,10)
#Look at the plot
par(mfrow=c(2,3))

hist(gibb.sample[-(1:burnIn),2],nclass=30, main= "posterior of beta1", xlab="True value = red line")
abline(v=mean(chain[-(1:burnIn),2]))
abline(v=beta[2],col="red")

hist(gibb.sample[-(1:burnIn),3],nclass=30, main= "posterior of sigma.sq", xlab="True value = red line")
abline(v=mean(chain[-(1:burnIn),3]))
abline(v=sigma.sq,col="red")

hist(gibb.sample[-(1:burnIn),4],nclass=30, main= "posterior of nu.sq", xlab="True value = red line")
abline(v=mean(chain[-(1:burnIn),4]))
abline(v=nu.sq,col="red")

plot(ts(gibb.sample[-(1:burnIn),2]),main="beta1")
abline()
plot(ts(gibb.sample[-(1:burnIn),3]),main="sigma.sq")
plot(ts(gibb.sample[-(1:burnIn),4]),main="nu.sq")
