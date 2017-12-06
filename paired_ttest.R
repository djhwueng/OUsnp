rm(list=ls())
library(phytools)
library(ape)
# rbinom(n=100,size=1,prob = 0.5) 暫時用白努力分佈造出隨機得01資料來造樹
tree2<-"(1:0,((2:0.607275,(((3:0.450259,6:0.387621):7.003134,5:8.835091):6.509057,(4:7.430765,7:0.251266):0.233285):0.219600):7.353626,8:0.631373):0.653199);"
phy2<-read.newick(text = tree2)
plot(phy2)
V<-vcv(phy2)+1e-5
V

# 用例子分解V矩陣
# x<-c(4,12,-16)
# y<-c(12,37,-43)
# z<-c(-16,-43,98)
# k<-rbind(x,y,z)
# k
# U<-chol(k)
# solve(U)

# decoposition為k的特徵值
decom<-eigen(V)
decom$values


# PDP^t
decom$vectors%*%sqrt(diag(decom$values))%*%sqrt(diag(decom$values))%*%t(decom$vectors)
Q_inv<-decom$vectors%*%diag(1/sqrt(decom$values))
# 
# decom$vectors%*%diag(decom$values)%*%t(decom$vectors)
# decom$vectors%*%sqrt(diag(decom$values))%*%sqrt(diag(decom$values))%*%t(decom$vectors)
# Q<-decom$vectors%*%sqrt(diag(decom$values))
# Q_inv<-sqrt(diag(decom$values))%*%t(decom$vectors)
# Q%*%Q_inv
# 


#Z<-Q_inv%*%matrix(test.data,ncol=1)

  
# 檢查兩個不同分佈的雙母體t檢定是否顯著
s1<-rnorm(4, mean=5,sd=2)
s2<-rnorm(4,mean=3,sd=3)
s=c(s1,s2)
s
paired.t<-t.test(s1,s2,paired = TRUE)
ls(paired.t)
paired.t$p.value


test.data<-sample(c(s1,s2))
test.data
snp.site<-c(1,0,1,0,1,0,0,1)
s1<-test.data[which(snp.site==1)]
s2<-test.data[which(snp.site==0)]
result.t<-t.test(s1,s2,paired = TRUE)
result.t$p.value



# 做8個種族，每個種族有10000個snp的檢定
sim.data<-matrix(rbinom(80000,size = 1,prob = 0.5), nrow = 8)
# ?rbinom
head(sim.data)
sim.data[,1:10]
row.names(sim.data)<-c(LETTERS[seq(from=1,to=8)])
?letters
sim.data[,1:10]
dim(sim.data)


# 寫一個函數讓模擬跑出10000次的檢定並且印出p-value小餘0.05的位置
paired.t.snp<-function(snp=snp,trait=trait){
  snp<-snp.site
  trait<-test.data
  s1<-trait[which(snp==1)]
  s2<-trait[which(snp==0)]
  length(s1)
  length(s2)
  result.t<-t.test(s1,s2,paired = FALSE)
  return(result.t)
}


p.raw.array<-array(0,dim(sim.data)[2])
p.uni.array<-array(0,dim(sim.data)[2])
for(snpindex in 1:dim(sim.data)[2]){
  #snpindex
  snp.site<-sim.data[,snpindex]
  #snp.site
  #test.data
#  print(try(paired.t.snp(snp=snp.site,trait=test.data)))
  p.raw.array[snpindex]<-try(paired.t.snp(snp=snp.site,trait=test.data)$p.value)
  p.uni.array[snpindex]<-try(paired.t.snp(snp=snp.site,trait=Q_inv%*%matrix(test.data,ncol=1))$p.value)
  }

length(which(p.raw.array<0.05))
length(which(p.uni.array<0.05))


# which(p.raw.array<0.05)
#which(p.array<0.05)

