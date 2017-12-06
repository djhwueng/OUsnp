rm(list=ls())
library(ape)
#https://github.com/bbanbury/phrynomics
#library(phrynomics)
# setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/micezhangdataset")
setwd("~/Dropbox/YiWeiHsu/Rcode/main.code/micezhangdataset")
snpdata <- read.csv(file = "snps.dat.collapsed.csv", header = TRUE)
pheno <- read.csv(file = "pheno_raw.csv" , header = TRUE)
# pheno

chromosome1<-subset(snpdata, snpdata$chr=="chr1")
# chr <- snpdata[,1]
# chr1<-chr[chr=="chr1"]
# snpdata[chr1,]


#write.csv(chromosome1,file = "chromosome1.csv")
#chromosome1.1<-chromosome1[,-1:-4]
#dim(chromosome1.1)
#head(chromosome1.1)

#chr1<-t(chromosome1.1)
#dim(chr1)
#chr1[1:5,1:5]
#system("pwd")
#snpdata <-ReadSNP(chr1)
#snpdata
#cbind(chr1,2:4164)
#write.csv(chromosome1.1,file = "chromosome1.1")
#chr1.1<-read.csv(file = "chromosome1.1.csv", header = TRUE, sep = " ")
#dim(chr1)
#chr1<-t(chr1)
#chr1[1:5,1:5]

#############造chr1.1fasta##############
chr1.1<-as.matrix(read.csv(file = "chr1.1.csv" , header = FALSE))
chr1.1.array<-array(0, c(dim(chr1.1)[1],2))
chr1.1.array[,1]<-chr1.1[,1]
chr1.1.array[,1]
for(Index in 1:dim(chr1.1)[1]){  
   chr1.1.array[Index,2]<-  as.character(paste(chr1.1[Index, 2:dim(chr1.1)[2]],collapse=""))
  }

chr1.1.array[1,2]
colnames(chr1.1.array)<-c("chr","seq")
chr1.1.array<-as.data.frame(chr1.1.array)
chr1.1.array$chr


chr1.1.array$seq[2] 


chr1.1.array[1,]
system("pwd")
fa = character(2 * nrow(chr1.1.array)) 
fa[c(TRUE, FALSE)] = sprintf("> %s", chr1.1.array$chr) 
fa[c(FALSE, TRUE)] = as.character(chr1.1.array$seq)
fa
writeLines(fa, "chr1.1.fasta") 

######################################

chr1.0<-as.matrix(read.csv(file = "chr1.0.csv" , header = FALSE))
chr1.0.array<-array(0, c(dim(chr1.0)[1],2))
chr1.0.array[,1]<-chr1.0[,1]
chr1.0.array[,1]
for(Index in 1:dim(chr1.0)[1]){  
  chr1.0.array[Index,2]<-  as.character(paste(chr1.0[Index, 2:dim(chr1.0)[2]],collapse=""))
}

chr1.0.array[1,2]
colnames(chr1.0.array)<-c("chr","seq")
chr1.0.array<-as.data.frame(chr1.0.array)
chr1.0.array$chr


chr1.0.array$seq[2] 


chr1.0.array[1,]
system("pwd")
fa = character(2 * nrow(chr1.0.array)) 
fa[c(TRUE, FALSE)] = sprintf("> %s", chr1.0.array$chr) 
fa[c(FALSE, TRUE)] = as.character(chr1.0.array$seq)
fa
writeLines(fa, "chr1.0.fasta") 




##############################################

# library(seqinr)
# chr1.1fas<-read.fasta("chr1.1.fasta")
# library(ape)
# write.nexus.data(chr1.1fas, file="chr1.1.nex" )


# ?write.nexus.data

#Use this website to get a string for the snp
#https://stackoverflow.com/questions/2098368/concatenate-a-vector-of-strings-character

# sdata = c('a','b','c')
# paste(sdata[1],sdata[2],sdata[3],sep='')
# ??pasta
# paste(sdata,collaose='')



chromosome5<-subset(snpdata, snpdata$chr=="chr5")
write.csv(chromosome5,file = "chromosome5.csv")


chr5.0<-as.matrix(read.csv(file = "chr5.0.csv" , header = FALSE))
chr5.0<-t(chr5.0)
chr5.0.array<-array(0, c(dim(chr5.0)[1],2))
chr5.0.array[,1]<-chr5.0[,1]
chr5.0.array[,1]
for(Index in 1:dim(chr5.0)[1]){  
  chr5.0.array[Index,2]<-  as.character(paste(chr5.0[Index, 2:dim(chr5.0)[2]],collapse=""))
}

chr5.0.array[1,2]
colnames(chr5.0.array)<-c("chr","seq")
chr5.0.array<-as.data.frame(chr5.0.array)
chr5.0.array$chr


chr5.0.array$seq[2] 


chr5.0.array[1,]
system("pwd")
fa = character(2 * nrow(chr5.0.array)) 
fa[c(TRUE, FALSE)] = sprintf("> %s", chr5.0.array$chr) 
fa[c(FALSE, TRUE)] = as.character(chr5.0.array$seq)
fa
writeLines(fa, "chr5.0.fasta") 

######################################################################
chr5.1<-as.matrix(read.csv(file = "chr5.1.csv" , header = FALSE))
chr5.1<-t(chr5.1)
chr5.1.array<-array(0, c(dim(chr5.1)[1],2))
chr5.1.array[,1]<-chr5.1[,1]
chr5.1.array[,1]
for(Index in 1:dim(chr5.1)[1]){  
  chr5.1.array[Index,2]<-  as.character(paste(chr5.1[Index, 2:dim(chr5.1)[2]],collapse=""))
}

chr5.1.array[1,2]
colnames(chr5.1.array)<-c("chr","seq")
chr5.1.array<-as.data.frame(chr5.1.array)
chr5.1.array$chr

chr5.1.array$seq[2] 

chr5.1.array[1,]
system("pwd")
fa = character(2 * nrow(chr5.1.array)) 
fa[c(TRUE, FALSE)] = sprintf("> %s", chr5.1.array$chr) 
fa[c(FALSE, TRUE)] = as.character(chr5.1.array$seq)
fa
writeLines(fa, "chr5.1.fasta") 

############################################################################
chromosome8<-subset(snpdata, snpdata$chr=="chr8")
write.csv(chromosome8,file = "chromosome8.csv")


chr8.1<-as.matrix(read.csv(file = "chr8.1.csv" , header = FALSE))
chr8.1<-t(chr8.1)
chr8.1.array<-array(0, c(dim(chr8.1)[1],2))
chr8.1.array[,1]<-chr8.1[,1]
chr8.1.array[,1]
for(Index in 1:dim(chr8.1)[1]){  
  chr8.1.array[Index,2]<-  as.character(paste(chr8.1[Index, 2:dim(chr8.1)[2]],collapse=""))
}

chr8.1.array[1,2]
colnames(chr8.1.array)<-c("chr","seq")
chr8.1.array<-as.data.frame(chr8.1.array)
chr8.1.array$chr

chr8.1.array$seq[2] 

chr8.1.array[1,]
system("pwd")
fa = character(2 * nrow(chr8.1.array)) 
fa[c(TRUE, FALSE)] = sprintf("> %s", chr8.1.array$chr) 
fa[c(FALSE, TRUE)] = as.character(chr8.1.array$seq)
fa
writeLines(fa, "chr8.1.fasta") 
#########################################################################################
chr8.0<-as.matrix(read.csv(file = "chr8.0.csv" , header = FALSE))
chr8.0<-t(chr8.0)
chr8.0.array<-array(0, c(dim(chr8.0)[1],2))
chr8.0.array[,1]<-chr8.0[,1]
chr8.0.array[,1]
for(Index in 1:dim(chr8.0)[1]){  
  chr8.0.array[Index,2]<-  as.character(paste(chr8.0[Index, 2:dim(chr8.0)[2]],collapse=""))
}

chr8.0.array[1,2]
colnames(chr8.0.array)<-c("chr","seq")
chr8.0.array<-as.data.frame(chr8.0.array)
chr8.0.array$chr

chr8.0.array$seq[2] 

chr8.0.array[1,]
system("pwd")
fa = character(2 * nrow(chr8.0.array)) 
fa[c(TRUE, FALSE)] = sprintf("> %s", chr8.0.array$chr) 
fa[c(FALSE, TRUE)] = as.character(chr8.0.array$seq)
fa
writeLines(fa, "chr8.0.fasta") 
