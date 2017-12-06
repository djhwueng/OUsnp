#install.packages("phyclust")
#install.packages("phangorn")
#install.packages("phytools")
# library(phytools)
# library(phyclust)
# library(phangorn)
#demo("toy",package="phyclust")

#vcv.phylo
# 
# X<-seq.data.toy$org
# X
# X.class<-as.numeric(gsub(".*-(.)", "\\1",seq.data.toy$seqname))
# dim(X)
# 
# #windows()
# plotdots(X,X.class)
# X
# ret <- phyclust.edist(X,edist.model = .edist.model[3])
# length(ret)
# ?phyclust.edist
# 
# upgma_tree<-upgma(ret)
# upgma_tree
# plot(upgma(ret))
# 
# 

#https://en.wikipedia.org/wiki/UPGMA
#ret.tree<-nj(ret)
#ret.tree$tip.label
#ls(ret.tree)
#plot(root(ret.tree,outgroup="50"))
#dim(X)
#X[1:10,1:10]
#plotnj(ret.tree,X.class=X.class)
#?plotnj

#EMC.2 <- .EMControl(init.procedure="emEM")
#EMC.2
#set.seed(1234)
#(ret.2<-phyclust(X,4,EMC=EMC.2))
#RRand(ret.2$class.id,X.class)
#.code.type
#.nucleotide
#.snp
#?ms

# Y <- rnorm(100 )
# library(ape)
# V<-vcv.phylo(upgma_tree)
# ls(upgma_tree)
# upgma_tree$edge
# dim(V)
# 
# max(V)
# V[1:5,1:5]/max(V)
# Va(alpha=0.001, V=V)[1:5,1:5]
# 
# 
# 
# #DNA seq
# X
# #trait data
# Y
#######################################################
###############   CODE FOR SNP, paird t-test ############### 
#######################################################
Y1<-rnorm(50,mean=6)
Y2<-rnorm(50, mean=3)
Y<-sample( c(Y1,Y2))
D<-matrix(rbinom(1000*100,1,prob=0.5),ncol=1000,nrow=100)
head(D)
D[1:100,1:10]
D[,1]
p.value.array<-array(1,c(1000))

for(index in 1:1000){
    nY1<-Y[D[,index]==0]
    nY2<-Y[D[,index]==1]
    p.value.array[index]<-t.test(nY1,nY2)$p.value
    }
#######################################################
###############   CODE FOR SNP, D, VCV  ############### 
#######################################################
rm(list=ls())
library(phytools)
library(phyclust)
library(phangorn)

  Va<-function(alpha,V=V){
      V<-V/max(V)
      A<-exp(-2*alpha*(1-V))
      B<- (1-exp(-2*alpha*V))/(2*alpha)  
      return(A*B)
      }
  #Finish example
  phy<- read.newick(text="(((F:0.06,E:0.06)G:0.28,D:0.34)I:0.07,(C:0.35,(A:0.12,B:0.12)H:0.23)J:0.06)k:0;")
  plot(phy)
  #http://www.phytools.org/eqg/Exercise_3.2/
  nodelabels()
  tiplabels()
  edgelabels(phy$edge.length)
  #At each SNP site, the tree can be partitioned into k clusters 
  #using only the earliest k-1 edges in the tree. 
  #cluster_tree<-function(k,tree=tree){
  #for n tips, a bifurcated tree has 
  #tree<-tree.newick
  ntip<-length(phy$tip.label)
  phy<-reorder(phy,"pruningwise")
  anc<-phy$edge[,1]
  des<-phy$edge[,2]
  N<-dim(phy$edge)[1]
  x <- rep(0, ntip + phy$Nnode)
  ROOT<-ntip+1
  x[ROOT]<-1
  relation<-cbind(phy$edge,matrix(phy$edge.length,ncol=1))
  colnames(relation)<-c("anc","des","length")

  k <- 2  #number of clusters
  D<-array(0,c(ntip,k))
  D
  #find the two ancestors firsts which are 8 and 10
  #find descendants and assign them to each cluster separately
  Inode <- subset(des,des>ntip)
  Inode
  count_cluster<-0
  for(i in length(Inode): (length(Inode)-k +1)){
#  assign( paste("C",i,seq=""), c(Inode[i],getDescendants(phy, node=Inode[i])))
   desc.temp<-getDescendants(phy, node=Inode[i])
        assign(paste("I",Inode[i],sep=""), subset(desc.temp,desc.temp<=ntip))
    #  C2<-getDescendants(phy, node=des[N-1])
   print(get(paste("I",Inode[i],sep="")))
   count_cluster<-count_cluster+1 
  }
  
  I9
  I10
  I11
  I8
count_cluster
  for(clusterIndex_outer in length(Inode): (length(Inode) -count_cluster +2)){
 #   clusterIndex_outer<-4
    #clusterIndex_outer
    Io<-get(paste("I",Inode[clusterIndex_outer],sep=""))
    #Io
      for(clusterIndex_Inner in (clusterIndex_outer-1):(length(Inode) -count_cluster+1)){
#        clusterIndex_Inner<-(length(Inode) -count_cluster +1)
#        clusterIndex_Inner
        In<-get(paste("I",Inode[clusterIndex_Inner],sep=""))
    #    if(clusterIndex_Inner != clusterIndex_outer){
#    In
          if (all(In %in% Io)){
          assign(paste("I",Inode[clusterIndex_outer],sep="") ,setdiff(Io,In))
#          }
        }    
      }
    }
  
  I9
  I10
  I11
  I8
  
  
  
  Inode_count<-0
  tree.D<-array(0,c(ntip,k ))
  for(i in length(Inode): (length(Inode)-count_cluster+1) ){
  
    print(get(paste("I",Inode[ i ],sep="")))
    d.index <-get(paste("I",Inode[i],sep=""))
    d.index
    Inode_count<-Inode_count+1
    D[d.index,Inode_count]<-1
    #tree.D[d.index,Inode_count]<- phy$edge.length[which(phy$edge[,2]==Inode[i])]
    tree.D[d.index,Inode_count]<- nodeHeights(phy)[,2][which(phy$edge[,2]==Inode[i])]
    #i<-4
    }
  D
  colnames(D)<- Inode[length(Inode): (length(Inode)-count_cluster+1)  ]
  D
  
  
  nodeHeights(phy)
  des
  phy$edge
  tree.D
  
  t(tree.D)
  
  V1<-sqrt(tree.D%*%t(tree.D))
  #unique(tree.D)
  #?unique
  #THINK ABOUT HOW TO MAKE THE RIGHT.EXACT MATRIX HERE.
  V1
   
 
  
#  row.names(V1)<- 
  
  #Inode
  #check from root
  addheight1<-nodeHeights(phy)[,2][which(phy$edge[,1]==ntip+1)[1]]
  addheight1
  D_1<- getDescendants(phy, phy$edge[which(phy$edge[,1]==ntip+1),2][1])
  D_1<-subset( D_1,D_1<=ntip    )
  D_1
  
  putD1<-D_1[which(apply(V1[D_1,]!=0,1,sum)==1)]  
  V1[putD1,D_1]<-addheight1
  V1[D_1,putD1]<-addheight1
  V1

  
  addheight2<-nodeHeights(phy)[,2][which(phy$edge[,1]==ntip+1)[2]]
  addheight2
  D_2<- getDescendants(phy, phy$edge[which(phy$edge[,1]==ntip+1),2][2])
  D_2<-subset( D_2,D_2<=ntip    )
  D_2
  
  putD2<-D_2[which(apply(V1[D_2,]!=0,1,sum)==1)]  
  V1[putD2,D_2]<-addheight2
  V1[D_2,putD2]<-addheight2
  V1
  diag(V1)<-diag(vcv.phylo(phy))[1]
  print(V1)
  #GOOD CODE fOR TREE VCV, SO WE CAN MOVE fORWARD TO PUT PROCESS ON TREE FOR TRAIT EVOLUTION 
  ########################################################################
  # we have D and V now, so we want to do parameter estimation on the mu, sigma, and alpha
  # the first thing is to do brownian motion
  # so given trait, tree, and k, we will report the estimate of mu and sigma first by the formula (4) in 
  #(4) and (5) in thompson

  
  
  # phy$edge
  # 
  # 
  # 
  # colnames(D)<-Inode
  # D
  # #CURRENT CODE WORKS FOR AN EXAMPLE
  # # we now have a correct structure for D, we next want to
  # 
  # 
  # phy$edge.length[which(phy$edge[,2]==Inode[4])]
  # 
  # 
  # D
  
  
  
  
  # 
  # #Next we want to add tip and then drop tip
  # #first look at what is the tree
  # new.tree<- phy
  # #http://blog.phytools.org/2012/11/adding-single-tip-to-tree.html
  # 
  # #http://blog.phytools.org/2013/01/adding-new-tips-at-random-to-phylogeny.html
  # plot(new.tree)
  # edgelabels(phy$edge.length)
  # nodelabels()
  # tiplabels()
  
  # new.tree$edge
  # relation
  # 
  # new.tree<- phy#pbtree(n=6)
  # 
  # new.tree$edge
  # phy$edge
  # 
  # new.tree$edge.length[which(des==3)]
  # 
  # new.tree$edge[,2]
  # des
  # 
  # des
  # position<-new.tree$edge.length[which(tree$edge[,2]==3)]
  
  
  
  # new.tree$edge.length[which(tree$edge[,2]==3)]
  # 
  # plot(bind.tip(new.tree,"t01", where=5, position=0.1))
  # 
  # new.tree2<-bind.tip(new.tree,"t01", where=3, position=0.339999)
  # 
  # plot(new.tree2)
  # edgelabels(phy$edge.length)
  # new.tree$edge.length
  # 
  # ls(new.tree)
  # 
  
  
  # 
  # tree <- pbtree(n=6)
  # plot(tree)
  # tiplabels()
  # nodelabels()
  # 
  # d1.tree<-drop.tip(tree,1)
  # plot(d1.tree)
  # nodelabels()
  # tiplabels()
  # node<-sample(c(1:length(tree$tip),2:tree$Nnode+length(tree$tip)),size = 1)
  # node<-9
  # position <- tree$edge.length[which(tree$edge[,2]==node)]
  # new.tree<-bind.tip(tree,"t02",where=node,position = position)
  # plot(new.tree)
  # 
  # 
  # 
  # class(tip)<-"phylo"
  # btree<-bind.tree(tree,tip,where=16)
  # plotTree(btree , node.numbers = T)
  # 
  # 
  # 
  # 
  # relation
  # ?assign
  #   
  # tips( phy,node=7)
  # 
  # splitTree(phy,split=7)
  # 
  # ?splitTree
  # 
  # ?getDescendants
  # subset(C1,C1<=ntip)
  # subset(C2,C2<=ntip)
  # D<-array(0,c(ntip,k))
  # D[subset(C1,C1<=ntip),1]<-1
  # D[subset(C2,C2<=ntip),2]<-1
  # D
  # 
  
  
  
  # for(i in N:1){
  # #  i<-N-9
  # #  print(des[i])
  # #  print(anc[i])
  #   x[des[i]]<- x[anc[i]]+ rnorm(1,sd=phy$edge.length[i])
  # 
  #   
  #   
  #   #    x
  #       #print( c(anc[i],des[i] ))
  #   #print(  x[des[i]])
  #   }
  # 
#   
#   
# #  http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
#  tree <- pbtree(n=10)
#  tree$tip.label<-1:10
#  par(mar=c(1,1,1,1))
#  plot(tree,font = 1)
#  nodelabels()
#   getDescendants(tree, node = 17)
#   getDescendants(tree, node = 12)
#   getDescendants(tree, node = 19)
#   x
  
#   if(k==1){newtree=tree
#     D<-array(1,c(n,k))
#     }
#   #remove tips and then add tips
#   
# #    return(list(D=D,newtree=newtree))
# #  }
# 
# 
# V1<-vcv.phylo(tree.newick)
# tree.newick.1<- reorder(tree.newick , "pruningwise")
# plot(tree.newick.1)
# tree.newick.1$edge
# 
# plot(extract.clade(tree.newick.1,8))
# 
# nodelabels()
# tiplabels()
# 
# 
# data(bird.families)
# 
# tree.newick.1$edge
# tree.newick.1$edge.length
# 
# tt10<-extract.clade(tree.newick,10)
# plotTree(tt10)
# 
# ?reorder
# length(tree.newick$tip.label)
# plot(tree.newick)
# matrix(tree.newick$edge.length,ncol=1)
# cbind(tree.newick$edge,matrix(tree.newick$edge.length,ncol=1))
# str(tree.newick)
# 
