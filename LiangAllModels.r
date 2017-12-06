numModels<-16
 
TableMatrix<-function(table){
    out.table<-array(0,c(dim(table)))
    for(i in 1:dim(table)[1]){
        for (j in 1:dim(table)[2]){
            out.table[i,j]<-table[i,j]
          }
        }
    return(out.table)
    }
       
doSimAnalysisReplicate<-function(treetype=treetype,datatype=datatype,ntax=ntax,nsite=nsite,outer.Index=outer.Index){
    nameString<-paste(treetype,datatype,ntax,nsite,outer.Index,sep="_")
    newdir<-paste("/Users/Shared/Brian.Tony/Brian.Liang.first_run/",nameString,sep="")
    system(paste("mkdir ",newdir))
    setwd(newdir)
    
    phy<-compute.brlen(stree(ntax,type=treetype))
   # http://bioweb2.pasteur.fr/docs/seq-gen/

    if(datatypeIndex==1){#simulate with JC  
       setup<-paste("-mJC -l",nsite," -t1 -s0.1 -d1 -f0.25,0.25,0.25,0.25 -on", sep="")
       data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }

    if(datatypeIndex==2){#simulate with HKY     
       setup<-paste("-mHKY -l",nsite," -t2 -s0.1 -d1 -f0.32,0.30,0.10,0.27 -on", sep="")
       data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }


    if(datatypeIndex==3){#simulate with GTR 
     setup<-paste("-mGTR -l",nsite," -r15.87,33.84,7.36,2.13,43.32,1 -s0.1-d1 -f0.32,0.30,0.11,0.28 -on", sep="")
      data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }

    if(datatypeIndex==4){#simulate with HKY + gamma  rate    
       setup<-paste("-mHKY -l",nsite," -t2 -s0.1 -d1 -a0.5 -f0.32,0.30,0.10,0.27 -on", sep="")
       data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }

    if(datatypeIndex==5){#simulate with GTR + gamma + I rate 
       setup<-paste("-mGTR -l",nsite," -r15.87,33.84,7.36,2.13,43.32,1 -s0.1-d1 -f0.32,0.30,0.11,0.28 -i0.3 -a0.5 -on", sep="")
       data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }
   
    system(paste("cp /Users/Shared/Brian.Tony/Brian.Liang/modelblockKLprerun.nex ."))
    system("paup -n modelblockKLprerun.nex > /dev/null ")

    system(paste("cp /Users/Shared/Brian.Tony/Brian.Liang/parse_paupOutKL.pl ."))
    system ("perl parse_paupOutKL.pl")
    
    outer.table<- TableMatrix((read.table("summary.txt",stringsAsFactors=FALSE,header=TRUE)[1:numModels,2:13]))
    for(outer.tableIndex in 1: length(outer.table[,1])){
       if(outer.table[,1][outer.tableIndex]=="infinity"){outer.table[,1][outer.tableIndex]<-"Inf"}   
       }
       
    new.outer.table<-array(0,c(dim(outer.table)))
    for(i in 1:dim(outer.table)[1]){
      for(j in 1:dim(outer.table)[2]){
        new.outer.table[i,j]<-as.numeric(outer.table[i,j])
        }
      }
    
    outer.table<-new.outer.table
    
    outer.table[which(outer.table==0)]<-0.0000000001

   
    system("perl innerKL.pl > /dev/null ")
    system("rm model.scores")
    #Now we are going to inner loop, you will use the Q and bf for optimization.
    for(theta.Index in 1:innerloop){
       #print(c(datatypeVect[datatypeIndex] ,nsite,outer.Index,theta.Index))
       
       if(datatypeIndex==1){#simulate with JC  
       setup<-paste("-mJC -l",nsite," -t1 -s0.1 -d1 -f0.25,0.25,0.25,0.25 -on", sep="")
       data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }

    if(datatypeIndex==2){#simulate with HKY     
       setup<-paste("-mHKY -l",nsite," -t2 -s0.1 -d1 -f0.32,0.30,0.10,0.27 -on", sep="")
       data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }


    if(datatypeIndex==3){#simulate with GTR 
     setup<-paste("-mGTR -l",nsite," -r15.87,33.84,7.36,2.13,43.32,1 -s0.1-d1 -f0.32,0.30,0.11,0.28 -on", sep="")
      data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }

    if(datatypeIndex==4){#simulate with HKY + gamma  rate    
       setup<-paste("-mHKY -l",nsite," -t2 -s0.1 -d1 -a0.5 -f0.32,0.30,0.10,0.27 -on", sep="")
       data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }
    if(datatypeIndex==5){#simulate with GTR + gamma + I rate 
       setup<-paste("-mGTR -l",nsite," -r15.87,33.84,7.36,2.13,43.32,1 -s0.1-d1 -f0.32,0.30,0.11,0.28 -i0.3 -a0.5 -on", sep="")
       data<-seqgen(opts = setup,rooted.tree = phy)
       write(data,"data.nex")
       }
    
     system("paup -n innerKL.nex > /dev/null ")

     system ("perl parse_paupOutKL.pl")

      inner.table<- TableMatrix((read.table("summary.txt",stringsAsFactors=FALSE,header=TRUE)[1:numModels,2:13]))
      
      for(inner.tableIndex in 1: length(inner.table[,1])){
       if(inner.table[,1][inner.tableIndex]=="infinity"){inner.table[,1][inner.tableIndex]<-"Inf"}   
       }
      
      new.inner.table<-array(0,c(dim(inner.table)))
      for(i in 1:dim(inner.table)[1]){
      for(j in 1:dim(inner.table)[2]){
        new.inner.table[i,j]<-as.numeric(inner.table[i,j])
        }
      }
      inner.table<-new.inner.table
      
      inner.table[,2:dim(inner.table)[2]]<-outer.table[,2:dim(outer.table)[2]]
      inner.table[which(inner.table==0)]<-0.0000000001
      #print(inner.table)
	  #print("before KL.mean.theta")
	 #print(KL.mean.theta)
       KL.mean.theta[,treetypeIndex,datatypeIndex,ntaxIndex,nsiteIndex]<-KL.mean.theta[,treetypeIndex,datatypeIndex,ntaxIndex,nsiteIndex]+inner.table[,1]
	#print("after KL.mean.theta")
       }#end of innerloop

    KL.mean.f[,treetypeIndex,datatypeIndex,ntaxIndex,nsiteIndex]<-KL.mean.theta[,treetypeIndex,datatypeIndex,ntaxIndex,nsiteIndex]/innerloop

    print(paste("Now saving: ",treetype,datatype,ntax,nsite,outer.Index,".Rsave"))
    try(save(KL.mean.f,AICarray,outer.table,file=paste(treetype,datatype,ntax,nsite,outer.Index,".Rsave",sep=""),compress=TRUE))
    setwd("/Users/Shared/Brian.Tony/Brian.Liang.first_run/")
    return(TRUE)
    }

###simulate sequences
library(ape)
library(rgl)
library(phyclust)
library(phangorn)
library(doMC)
library(foreach)
registerDoMC(100)

TransitionRateJC<- matrix(c(1,1,1,1,1,1),nrow=1)  
TransitionRateHKY<-matrix(c(1,2,1,1,2,1),nrow=1)
TransitionRateGTR<-matrix(c(15.87,33.84,7.36,2.13,43.32,1),nrow=1)
TransitionRateArray<-rbind(TransitionRateJC,TransitionRateHKY,TransitionRateGTR,TransitionRateHKY,TransitionRateGTR)
#TransitionRateArray<-TransitionRateGTR

BasefreqJC<-matrix(c(0.25,0.25,0.25,0.25),nrow=1)
BasefreqHKY<-matrix(c(0.27,0.23,0.18,0.32),nrow=1)
BasefreqGTR<-matrix(c(0.32,0.30,0.11,0.28) ,nrow=1 )
BasefreqArray<-rbind(BasefreqJC,BasefreqHKY,BasefreqGTR,BasefreqHKY,BasefreqGTR)
#BasefreqArray<-BasefreqGTR

treetypeVect=c("left","balanced")
datatypeVect<-c("JC","HKY","GTR","HKY+G","GTR+G+I")
ntaxVect=c(4)
nsiteVect<- c(1000)#,5000,10000,50000)

outerloop<-2
innerloop<-2

KL.mean.theta<-array(0,c(numModels,length(treetypeVect),length(datatypeVect),length(ntaxVect),length(nsiteVect) ) )
KL.mean.f<-array(0,c(numModels,length(treetypeVect),length(datatypeVect),length(ntaxVect),length(nsiteVect)) )

AICarray<-array(0,c(numModels,length(treetypeVect),length(datatypeVect),length(ntaxVect),length(nsiteVect)  ) )


for(treetypeIndex in 1:length(treetypeVect)){
   treetype<-treetypeVect[treetypeIndex]
  for(datatypeIndex in 1:length(datatypeVect)){
   datatype<-datatypeVect[datatypeIndex]
   transitionrate<-TransitionRateArray[datatypeIndex,]
   basefreq<-BasefreqArray[datatypeIndex,]
   for(ntaxIndex in 1:length(ntaxVect)){
      ntax<-ntaxVect[ntaxIndex] 
      for(nsiteIndex in 1:length(nsiteVect)){
         nsite<-nsiteVect[nsiteIndex] 
 
    foreach(outer.Index=c(1:outerloop)) %dopar% {doSimAnalysisReplicate(treetype,datatype,ntax,nsite,outer.Index)} 
 
   	}#end of siteIndex
     }#ntaxIndex 
   }#end of datatypeIndex
}#end of treetypeIndex




