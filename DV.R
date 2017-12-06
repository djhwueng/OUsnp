#rm(list=ls())
library(phytools)
library(phyclust)
library(ape)
library(phangorn)
library(corpcor)
library(phylobase)

DV<- function(k=k, tree=tree){
#    phy<-tree
#    rm(list=ls())
#    phy<-rtree(8)
    k<-2
    phy<-reorder(tree, "postorder") #reorder(x, order = c("preorder", "postorder","caldewise","pruningwise"))
#    k<-4
#    plot(phy)
#    nodelabels()
#    tiplabels()
    ntip<-length(phy$tip.label)

    #this is for tree with one outgroup
    Inode<- matrix(unique(c(phy$edge)[!(c(phy$edge) %in% (1:(ntip)))]),nrow=1)
    Inode.length<-matrix(unique(c(nodeHeights(phy))[!(c(phy$edge) %in% (1:(ntip)))]),nrow=1)
    Inode.length<-rbind(Inode,Inode.length)
    Inode.length <- Inode.length[,order(Inode.length[2,])]
    Inode.length #we have the info from the root tip to do cluster correctly
    geo.phy<-extractTree(phy) #from phylobase
    geo.phy
    for(rootupindex in 1:k){
      #rootupindex<-3
      #print(subset(geo.phy,node.subtree=Inode.length[1,rootupindex]))
      tip.num.set <- subset(geo.phy,node.subtree=Inode.length[1,rootupindex])@label
      #subtree.tip <- array(0,c(length(tip.num.set)))
      subtree.tip<-as.numeric(unlist(tip.num.set))
    #  for(sbt.index in 1:length(subtree.tip)){
#        subtree.tip[sbt.index]<-as.numeric(unlist(strsplit(tip.num.set[sbt.index],split="t"))[2]) #for rtree with label t1, t2, ...
#        }
        assign(paste("I", Inode.length[1,rootupindex],sep=""),subtree.tip)
#        print(paste("I",Inode.length[1,rootupindex],sep=""))
#        print(get(paste("I",Inode.length[1,rootupindex],sep="")))
      }#end for

    #this is for tree with two clades in root
    if(length(setdiff( get(paste("I",Inode.length[1,1],sep="")), get(paste("I",Inode.length[1,2],sep=""))))!=1){
      Inode<- matrix(unique(c(phy$edge)[!(c(phy$edge) %in% (1:(ntip+1)))]),nrow=1)
      Inode.length<-matrix(unique(c(nodeHeights(phy))[!(c(phy$edge) %in% (1:(ntip+1)))]),nrow=1)
      Inode.length<-rbind(Inode,Inode.length)
      Inode.length <- Inode.length[,order(Inode.length[2,])]
      Inode.length #we have the info from the root tip to do cluster correctly
      geo.phy<-extractTree(phy) #from phylobase
      for(rootupindex in 1:k){
        #print(subset(geo.phy,node.subtree=Inode.length[1,rootupindex]))
        tip.num.set <- subset(geo.phy,node.subtree=Inode.length[1,rootupindex])@label
        subtree.tip<-as.numeric(unlist(tip.num.set))
        #subtree.tip <- array(0,c(length(tip.num.set)))
        #for(sbt.index in 1:length(subtree.tip)){
        #  subtree.tip[sbt.index]<-as.numeric(unlist(strsplit(tip.num.set[sbt.index],split="t"))[2])
        #  }
          assign(paste("I", Inode.length[1,rootupindex],sep=""),subtree.tip)
#          print(paste("I",Inode.length[1,rootupindex],sep=""))
#          print(get(paste("I",Inode.length[1,rootupindex],sep="")))
        }#end for
      }#end if
    ##NEED TO FIX BELOW FIRST TO GET CORRECT CLUSTER
    # we check if the latter is the descedant of the former
    #assign to D from k to 1
    D<-array(0,c(ntip,k))
    #tree.D<-array(0,c(ntip,k))
    stack.pos<-c(1:ntip)
    for(tipdownindex in k:1){
      #tipdownindex<- 3
      des.pos<-get(paste("I",Inode.length[1,tipdownindex],sep=""))
#      des.pos
      if(length(intersect(stack.pos,des.pos))!=0){
      insert.pos<-intersect(stack.pos,des.pos)
      D[insert.pos,tipdownindex] <-1
      #for(insindex in 1:length(insert.pos)){
      #  tree.D[insert.pos[insindex],tipdownindex]<-nodeHeights(phy)[,2][which(phy$edge[,2]== insert.pos[insindex])]
      #  }
      }
      #print(insert.pos)

      stack.pos<-setdiff(stack.pos,des.pos)
      #stack.pos
      } #end of for loop
     colnames(D)<- Inode.length[1,1:k]
     D
     #colnames(tree.D)<-Inode.length[1,1:k]
     nonzero.column.index<-!apply(D,2,sum)==0
     #nonzero.column.index
     D<-D[,nonzero.column.index]

     #tree.D<-tree.D[,nonzero.column.index]
     #tree.D
#     print(D)
#     print(dim(D))


     #tree.D.vcv


     #drop.tip(phy,)
    # need to check apply(D,1,sum)==1
     #V<-tree.D%*%t(tree.D) # covariance for non-diagonal
     #diag(V)<-node.depth.edgelength(phy)[1:ntip] # need to take care the order with D
     #V<-sqrt(tree.D%*%t(tree.D)) #this is vcv for covariance matrix not including tip
     #vcv(phy)
     #node.depth.edgelength(phy)
#     return(list(D=D,V=V))
#     }
#     t(D)%*%vcv(phy)%*%D
     full.node.vcv<-vcvPhylo(phy)
     full.node.vcv
     #D
     #plot(phy)
     #nodelabels()
     #D
     #tree.D.anc.index<-

     ancIndex<- which(colnames(full.node.vcv) %in% names(nonzero.column.index))

     if(length(ancIndex) != k ){
      tree.D.anc.vcv<- array(0,c(k,k))
      tree.D.anc.vcv[1:length(ancIndex),1:length(ancIndex)] <- full.node.vcv[ancIndex,ancIndex]
      colnames(tree.D.anc.vcv)<-rownames(tree.D.anc.vcv)<- c(colnames(full.node.vcv)[ancIndex] ,setdiff(colnames(D),colnames(full.node.vcv)[ancIndex]))
    }else{tree.D.anc.vcv<-full.node.vcv[ancIndex,ancIndex]}

     #row.index
     #col.index
     D

     print(D)
     print(tree.D.anc.vcv)





     rownames(D)<-1:dim(D)[1]
     kmtx<-D
     for(kindex in 1:dim(D)[2]){
       #kindex<-4
       D
       kdiag <- diag(1,dim(tree.D.anc.vcv))
       rownames(kdiag)<-colnames(kdiag)<-rownames(tree.D.anc.vcv)
       kdiag
       delrow.index <- which(rownames(kdiag)== colnames(D)[kindex])
       delrow.index
       insrow<-array(0, c(sum(D[,kindex]), dim(tree.D.anc.vcv)[1]  ))
       colnames(insrow)<-colnames(kdiag)
       insrow[   ,delrow.index]<-1
       rownames(insrow)<- which(D[,kindex]==1)
       insrow
       if(delrow.index==1){ #only down
         kdown<-kdiag[(delrow.index+1):dim(kdiag)[1] ,  ]
         kdown
         kmtx<- rbind(insrow,kdown)
         kmtx
         tree.D.anc.vcv<-kmtx%*%tree.D.anc.vcv%*%t(kmtx)
         tree.D.anc.vcv
         }else if(delrow.index==dim(kmtx)[1]){#only up
         kup<- kdiag[1:(delrow.index-1),]
         kmtx<- rbind(kup,insrow)
         kmtx
         tree.D.anc.vcv<-kmtx%*%tree.D.anc.vcv%*%t(kmtx)
         tree.D.anc.vcv
         }else{
       D
       kup<- kdiag[1:(delrow.index-1),]
       kup
       if(is.null(dim(kup))){
         kup<-matrix(kup,nrow=1)
         rownames(kup)<-rownames(kdiag)[1]
         kup
         }
       kdown<-kdiag[(delrow.index+1):dim(kdiag)[1] ,  ]
       kdown
       if(is.null(dim(kdown))){
         kdown<-matrix(kdown,nrow=1)
         rownames(kdown)<-rownames(kdiag)[dim(kdiag)[1]]
         kdown
         }
       kmtx<- rbind(kup,insrow,kdown)
       kmtx
       tree.D.anc.vcv<-kmtx%*%tree.D.anc.vcv%*%t(kmtx)
       tree.D.anc.vcv
       D
       }#end of else
       #array(0,c(dim(tree.D.anc.vcv)[1]+ sum(D[,kindex])-1, dim(tree.D.anc.vcv)[2]))
       # colnames(D)[kindex]
       # tree.D.anc.vcv
       # rownames(Kmtx)<-which(D[,kindex]==1)
       # tree.D.anc.vcv<-Kmtx%*%tree.D.anc.vcv%*%t(Kmtx)
       # print(dim(tree.D.anc.vcv)[1])
       }#end of for kindex



     # tree.D.vcv<-D%*%tree.D.anc.vcv%*%t(D)
     diag(tree.D.anc.vcv)<- node.depth.edgelength(phy)[1:ntip]
     return(list(D=D,treevcv = tree.D.anc.vcv))
   }#end of function

      # D
     # print(tree.D.vcv)
     # print(vcv(phy))
     #D ok, tree not yet
     #tr<-rtree(10)
     #V<-vcv(tr)
     #vcv2phylo(V)

     #tree.D.vcv
     #vcv2phylo

     #tree.D
     #write tree first then add tip according to node and D.
     #finally return vcv.

     #http://blog.phytools.org/2013/01/




#http://blog.phytools.org/2012/11/adding-single-tip-to-tree.html
# tree2<-pbtree(n=10)
# plotTree(tree2,node.numbers=T)
#
# tree2$edge
# tip<-list(edge=matrix(c(2,1),1,2),tip.label="species name",edge.length=1.0,Nnode=1)
# class(tip)<-"phylo"
# tip
# btree<-bind.tree(tree2,tip,where=13)
# plot(btree)
# nodelabels()
#
# btree$edge
