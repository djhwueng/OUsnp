#rm(list=ls())
library(phytools)
library(phyclust)
library(ape)
library(phangorn)
library(corpcor)
library(phylobase)
library(TreeSim)
DV<- function(k=k, tree=tree){
    #tree <-  sim.bd.taxa.age(n=size, numbsim=1, lambda=1, mu=1, frac = 0.5, age=1, mrca = TRUE)[[1]]
    phy<-reorder(tree, "postorder") #reorder(x, order = c("preorder", "postorder","caldewise","pruningwise"))
    #plot(phy)
    #nodelabels()
    #tiplabels()
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
      tip.num.set <- subset(geo.phy,node.subtree=Inode.length[1,rootupindex])@label
      subtree.tip <- array(0,c(length(tip.num.set)))
      for(sbt.index in 1:length(subtree.tip)){
        subtree.tip[sbt.index]<-as.numeric(unlist(strsplit(tip.num.set[sbt.index],split="t"))[2]) #for rtree with label t1, t2, ...
        }
        assign(paste("I", Inode.length[1,rootupindex],sep=""),subtree.tip)
      }#end for
    #this is for tree with two clades in root
      # if(length(setdiff( get(paste("I",Inode.length[1,1],sep="")), get(paste("I",Inode.length[1,2],sep=""))))!=1){
      #  Inode<- matrix(unique(c(phy$edge)[!(c(phy$edge) %in% (1:(ntip+1)))]),nrow=1)
      #  Inode.length<-matrix(unique(c(nodeHeights(phy))[!(c(phy$edge) %in% (1:(ntip+1)))]),nrow=1)
      #  Inode.length<-rbind(Inode,Inode.length)
      #  Inode.length <- Inode.length[,order(Inode.length[2,])]
      #  Inode.length #we have the info from the root tip to do cluster correctly
      #  geo.phy<-extractTree(phy) #from phylobase
      #  for(rootupindex in 1:k){
      #    tip.num.set <- subset(geo.phy,node.subtree=Inode.length[1,rootupindex])@label
      #    subtree.tip <- array(0,c(length(tip.num.set)))
      #    for(sbt.index in 1:length(subtree.tip)){
      #      subtree.tip[sbt.index]<-as.numeric(unlist(strsplit(tip.num.set[sbt.index],split="t"))[2])
      #      }
      #      assign(paste("I", Inode.length[1,rootupindex],sep=""),subtree.tip)
      #    }#end for
      #  }#end if

    #assign to D from k to 1
    D<-array(0,c(ntip,k))
    stack.pos<-c(1:ntip)
    for(tipdownindex in k:1){
      des.pos<-get(paste("I",Inode.length[1,tipdownindex],sep=""))
     des.pos
      if(length(intersect(stack.pos,des.pos))!=0){
      insert.pos<-intersect(stack.pos,des.pos)
      D[insert.pos,tipdownindex] <-1
      }

      stack.pos<-setdiff(stack.pos,des.pos)
      } #end of for loop
     colnames(D)<- Inode.length[1,1:k]
     nonzero.column.index<-!apply(D,2,sum)==0
     #nonzero.column.index
     D<-D[,nonzero.column.index]
     full.node.vcv<-vcvPhylo(phy)
     ancIndex<- which(colnames(full.node.vcv) %in% names(nonzero.column.index))
     if(length(ancIndex) != k ){
      tree.D.anc.vcv<- array(0,c(k,k))
      tree.D.anc.vcv[1:length(ancIndex),1:length(ancIndex)] <- full.node.vcv[ancIndex,ancIndex]
      colnames(tree.D.anc.vcv)<-rownames(tree.D.anc.vcv)<- c(colnames(full.node.vcv)[ancIndex] ,setdiff(colnames(D),colnames(full.node.vcv)[ancIndex]))
    }else{tree.D.anc.vcv<-full.node.vcv[ancIndex,ancIndex]}
     rownames(D)<-1:dim(D)[1]
     kmtx<-D
     for(kindex in 1:dim(D)[2]){
       kdiag <- diag(1,dim(tree.D.anc.vcv))
       rownames(kdiag)<-colnames(kdiag)<-rownames(tree.D.anc.vcv)
       delrow.index <- which(rownames(kdiag)== colnames(D)[kindex])
       insrow<-array(0, c(sum(D[,kindex]), dim(tree.D.anc.vcv)[1]  ))
       colnames(insrow)<-colnames(kdiag)
       insrow[   ,delrow.index]<-1
       rownames(insrow)<- which(D[,kindex]==1)
       if(delrow.index==1){ #only down
         kdown<-kdiag[(delrow.index+1):dim(kdiag)[1] ,  ]
         if(is.null(dim(kdown))){
           kdown<-matrix(kdown,nrow=1)
           rownames(kdown)<-rownames(kdiag)[dim(kdiag)[1]]
         }
         kmtx<- rbind(insrow,kdown)
         tree.D.anc.vcv<-kmtx%*%tree.D.anc.vcv%*%t(kmtx)
         }else if(delrow.index==dim(kmtx)[1]){#only up
         kup<- kdiag[1:(delrow.index-1),]
         kmtx<- rbind(kup,insrow)
         tree.D.anc.vcv<-kmtx%*%tree.D.anc.vcv%*%t(kmtx)
         tree.D.anc.vcv
         }else{
       kup<- kdiag[1:(delrow.index-1),]
       if(is.null(dim(kup))){
         kup<-matrix(kup,nrow=1)
         rownames(kup)<-rownames(kdiag)[1]
         }
       kdown<-NULL
       if(delrow.index<dim(kdiag)[1]){
       kdown<-kdiag[(delrow.index+1):dim(kdiag)[1] ,  ]
       if(is.null(dim(kdown))){
         kdown<-matrix(kdown,nrow=1)
         rownames(kdown)<-rownames(kdiag)[dim(kdiag)[1]]
         }
       }

       kmtx<- rbind(kup,insrow,kdown)
       tree.D.anc.vcv<-kmtx%*%tree.D.anc.vcv%*%t(kmtx)
       }#end of else
       }#end of for kindex
     diag(tree.D.anc.vcv)<-  1 #node.depth.edgelength(phy)[1:ntip]

     return(list(D=D,V=tree.D.anc.vcv))
     }#end of function


# size<-20
# k=5
# tree <-  sim.bd.taxa.age(n=size, numbsim=1, lambda=1, mu=1, frac = 0.5, age=1, mrca = TRUE)[[1]]
# plot(tree)
# nodelabels()
# tiplabels()
# print(DV(k=k,tree=tree))
