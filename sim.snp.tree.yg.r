#given tree.size, sequence length, params.values
#return tree and trait
library(phytools)
library(gap)
library(phyclust)
library(OUwie)
library(ape)

#setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/ms/testms")
sim.snp.tree.yg<-function(taxa.size=taxa.size,seqlen=seqlen,ouparams=ouparams){
  sim.snp <- paste("ms ", taxa.size, " 1 -T -s ", seqlen, " > ms.out", sep="")
  system(sim.snp)
  msout1 <- read.ms.output("ms.out")
  snp<-matrix(unlist(msout1$gametes),nrow=taxa.size)
  write.nexus.data(snp,file="snp.nex",format="dna")
  #write.nexus.data(snp,file=paste("snp",seqlen,sep=""), format=dna)
  system("sed 's/DNA/standard/g' snp.nex >snp.1.nex")
  system("paup -n get.tree.snp.nex >/dev/null")
  system("grep 'PAUP_1' snp.tre | sed 's/^.*: //' > snp.1.tre ")
  temp.tree<-read.table(file="snp.1.tre")
  tree<-unlist(temp.tree[5])
  tree<-as.character(tree)
  write(tree,file="tree.txt")
  tree<-read.newick(file="tree.txt")

  alpha <- ouparams[1:2]
  sigma.sq <- ouparams[3:4]
  theta0 <- ouparams[5]
  theta <- ouparams[6:7]

  smp.tree <- sim.history(tree,Q=matrix(c(-1,1,1,-1),2,2))#Juke-Cantor
  sim.data <- OUwie.sim(smp.tree, simmap.tree=TRUE, scaleHeight=FALSE, root.age=0, alpha=alpha,sigma.sq=sigma.sq,theta0=theta0,theta=theta)
  return(list(tree=tree, trait=sim.data))
  }

#setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/BMsnp")

# taxa.size<-5
# seqlen<-100
# ouparams<-c(5,5,40,40,90,80,100)
# names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
#
# sim.tree.trait<- sim.snp.tree.yg(taxa.size=taxa.size,seqlen=seqlen,ouparams=ouparams)
# print(sim.tree.trait)
