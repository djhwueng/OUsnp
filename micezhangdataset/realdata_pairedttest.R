setwd("~/Dropbox/YiWeiHsu/Rcode/main.code/micezhangdataset")
load("contree1.RData")
load("contree5.RData")
load("contree8.RData")
load("realtree.RData")
chr1.0<-read.csv("chr1.0.csv" , header = T)
chr1.1<-read.csv("chr1.1.csv" , header = T)
chr5.0<-read.csv("chr5.0.csv" , header = T)
chr5.1<-read.csv("chr5.1.csv" , header = T)
chr8.0<-read.csv("chr8.0.csv" , header = T)
chr8.1<-read.csv("chr8.1.csv" , header = T)



s1<-which(chr1.0==1)
s2<-which(chr1.0==0)
result.t1.0<-t.test(s1,s2,paired = F)
ls(result.t1.0)
result.t1.0$statistic
s1

chr1.1<-read.csv("chr1.1.csv" , header = T)
s3<-which(chr1.1==1)
s4<-which(chr1.1==0)
result.t1.1<-t.test(s3,s4,paired = F)
result.t1.1$statistic
(result.t1.0$statistic+result.t1.1$statistic)/2
pheno<-read.csv("pheno_raw.csv" , header = T)
snps<- read.csv("snps.dat.collapsed.csv" , header = T)
#############################################################################
# sample(pheno$SBP)
dim(chr1.0)
N=4165
temp <-c()

for (i in 2:N){

  s1<-pheno$SBP[chr1.0[,i]==1]
  s2<-pheno$SBP[chr1.0[,i]==0]

  if(length(s1)>1){
    result.t1.0<-t.test(s1,s2,paired = F)
    temp[i]<-result.t1.0$statistic
    }else{temp[i]=0}
  }
temp
mean(temp[2:4165])
plot(temp[200:299])


N=4165
temp1.1<-c()

for (i in 2:N){

  s3<-pheno$SBP[chr1.1[,i]==1]
  s4<-pheno$SBP[chr1.1[,i]==0]

  if(length(s3)>1 & length(s4)>1){
    result.t1.1<-t.test(s3,s4,paired = F)
    temp1.1[i]<-result.t1.1$statistic
  }else{temp1.1[i]=0}
}
temp1.1
mean(temp1.1[2:4165])
plot(temp1.1[200:299])


(mean(temp[2:4165])+mean(temp1.1[2:4165]))/2




########################################################################################
i=2
s1<-pheno$SBP[chr1.0[,i]==1]
length(s1)
s2<-pheno$SBP[chr1.0[,2]==0]
length(s2)
result.t1.0<-t.test(s1,s2,paired = F)
result.t1.0$statistic
