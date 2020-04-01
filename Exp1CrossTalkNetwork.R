library(data.table)
library(plyr)
library(tidyverse)
#Learn where data files are stored.
source("dataConfig.R")

#Potentially, we can compute similarity of two pathways by using different gene similarity measures
simFiles  = c(bioProcessSimFile,molecularSimFile, proDomainSimFile, subcellularSimFile)

#cosineSim computes the standard cosine similarity of two gene vectors. Each entry in the vector is a pathway
#this assumes that vec1 and vec2 do not contain duplicate values.
cosineSim<-function(vec1,vec2){
   result = 0.0;
   if(is.na(vec1)||length(vec1)==0||is.na(vec2)||length(vec2)==0){
      result= 0.0;
   }
   
   inter = length(intersect(vec1,vec2))
   if(inter==0){ result =0.0}
   else{
      result = inter/(sqrt(length(vec1))*sqrt(length(vec2)))
   }
   return(result)
   
}

#GeneSim computes similarity of two genes by analyzing the pathways they appear in
# For a pair of genes, We use cosine similarity of their pathway vectors.
geneSim<-function(f=bioProcessSimFile){
   message("Similarity file is ",f)
   dat<- read.csv(f,header=F,sep=",",stringsAsFactors = F)
   colnames(dat)<-c("path","genelist")
   idData <- transform(dat,id=as.numeric(factor(path)))
   #how many genes are there?
   genes<-unique(unlist(as.list(strsplit(dat$genelist, "\t"))))
   numGenes = length(genes)
   #how many pathways?
   numPathways = nrow(dat)
   message(numGenes," genes, ",numPathways," pathways")
   newDat<-data.frame(gene=NA,pathway=NA)
   ps <- numeric()
   gs <- character()
   i<-0
   for(row in 1:nrow(idData)){
      p<-idData[row,"id"]
      for(g in unlist(as.list(strsplit(idData[row,"genelist"], "\t")))){
         i=i+1
         ps[i] <- p
         gs[i] <- g
      }
   }
   newDat<-data.frame(gs, ps, stringsAsFactors=FALSE)
   simMatrix = matrix(nrow=numGenes,ncol=numGenes)
   
   vs<-list()
   for(i in 1:numGenes){
      g1<-genes[i]
      vs[[g1]]<-newDat[newDat$gs==g1,]$ps
   }
   
   l1<-character()
   l2<-character()
   l3<-numeric()
   index=0
   seen<-list()
   for(row in 1:nrow(idData)){
      p<-idData[row,"id"]
      gList<-unlist(as.list(strsplit(idData[row,"genelist"], "\t")))
      for(i in 1:length(gList)){
         g1<-gList[i]
         vec1 = vs[[g1]]
         for(j in (i+1):length(gList)){
            
            g2<-gList[j]
               vec2 = vs[[g2]]
               index=index+1
               l1[index] <- g1
               l2[index] <- g2
               l3[index] <- cosineSim(vec1,vec2)
         } 
      }
   }
   simDat<-data.frame(l1,l2, l3, stringsAsFactors=FALSE)
   return(unique(simDat)) 
}



for (simFile in simFiles) {
   message(simFile," is being processed.")
   tryCatch({
      sim1<-readRDS(file = paste("results/",str_sub(simFile,start=-20),"geneSim.rds",sep=""))
      message("Previous similarity results are loaded")
   }, error = function(e) {
      message("Results were not found for the file. Gene similarity will be computed now.")
      result = geneSim(simFile)
      saveRDS(result,file = paste("results/",str_sub(simFile,start=-20),"geneSim.rds",sep=""))
   }, finally = {
      
   })
   
   
}
