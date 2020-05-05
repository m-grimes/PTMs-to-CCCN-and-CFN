library(data.table)
library(plyr)
library(tidyverse)
library(igraph)

#Learn where data files are stored.
source("dataConfig.R")

#Potentially, we can compute similarity of two pathways by using different gene similarity measures
simFiles  = c(bioProcessSimFile,molecularSimFile, proDomainSimFile, subcellularSimFile)

#cosineSim computes the standard cosine similarity of two gene vectors. 
#Each entry in the vector is a pathway
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
#
geneSim<-function(f=bioProcessSimFile){
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

computeEdgeWeights<-function(A,df){
   pList = unique(df$PATHWAY_ID)
   for (i in 1:length(pList)){
      for( j in (i+1):length(pList)){
         A[nrow(A) + 1,] = c(pList[[i]],pList[j],1)
      }
   }
   return(A)
}


for (simFile in simFiles) {
   message(simFile," is being processed to create gene similarity values.")
   tryCatch({
      sim1<-readRDS(file = paste("results/",str_sub(simFile,start=-20),"geneSim.rds",sep=""))
      message("Previously computed similarity results are loaded")
   }, error = function(e) {
      message("There is no previous save. Gene similarity will be computed now.")
      result = geneSim(simFile)
      saveRDS(result,file = paste("results/",str_sub(simFile,start=-20),"geneSim.rds",sep=""))
   }, finally = {
      
   })
}
   #rename the similarity columns
   colnames(sim1)<-c("gene1","gene2","similarity")
   #load pathways and their genes
   pathways = read.csv(pathwayFile)
   # how many genes are found in the pathways file?
   a<-(c(unique(sim1$gene1),unique(sim1$gene2)))
   b<-unique(pathways$GENE_SYMBOL)
   c<-intersect(a,b)
   message(length(c)," genes from ",pathwayFile," are also found in the similarity file ",simFile)
   
   #load clusters and their genes
   clusters = read.csv(clusterFile,sep="\t",stringsAsFactors = FALSE)
   #some rows contain multiple sites (e.g., TUBB4B ubi K379; TUBB2A ubi K379; TUBB2B ubi K379	1.225.171)
   #we will split those and add as rows
   s <- strsplit(clusters$Site, split = as.character(';'))
   clusters=data.frame( V1 = unlist(s),V2 = rep(clusters$Cluster, sapply(s, length)))
   colnames(clusters)<-c("Site","Cluster")
   #it is very easy to hate R
   clusters$Site<-trimws(as.character(clusters$Site))
   # we are interested in the symbol of the Site 
   # for example, we need to extract "TUBB4B" from TUBB4B ubi K379
   
   clusters$Site=sapply(strsplit(clusters$Site," "), `[`, 1)
   allPtmCounts = data.frame(table(clusters$Cluster))
   colnames(allPtmCounts)=c("Cluster","ptmCount")
   #now, we need to create pathways of clusters
   clusterPathways<- merge(clusters, pathways, by, by.x="Site", by.y="GENE_SYMBOL", sort = TRUE)
   clusterPathways$Cluster<-as.character(clusterPathways$Cluster)
   #how many ptms there are in each cluster
   invPtmCounts = data.frame(table(clusters$Cluster))
   colnames(invPtmCounts)=c("Cluster","ptmCount")
   
   #how many pathways are in each cluster
   invPwCounts = data.frame(table(clusterPathways$Cluster))
   colnames(invPwCounts)=c("Cluster","pathwayCount")
   
   #finally, we are ready for computations
   maxStep = 10
   percOfPathways = 2;
   dim = length(unique(clusterPathways$PATHWAY_ID))
   for (minPtmClusters in seq.int(50, 200, 10)){
      for (step in seq.int(6,maxStep,1)){
         
         percOfPtms = step / maxStep;
         graph = graph.empty( )
         A = data.frame( p1=character(), p2=character(),weight=numeric(), stringsAsFactors=FALSE)
         
         for(cl in unique(clusters$Cluster)){
            message(cl)
            message(proc.time())
            ptmCount = allPtmCounts[allPtmCounts$Cluster==cl,]$ptmCount
            proc.time()
            involvedPtmCount = invPtmCounts[invPtmCounts$Cluster==cl,]$ptmCount
            proc.time()
            involvedPathwayCount = invPwCounts[invPwCounts$Cluster==cl,]$pathwayCount
            proc.time()
            
            #rule 1
            if (ptmCount==0&&
                involvedPtmCount / (1.0 * ptmCount) < percOfPtms) {
               next;
            }
            #rule 2
            if (length(involvedPathwayCount)>0&&
                involvedPathwayCount / (involvedPtmCount) < percOfPathways) {
               next;
            }
            #left at here
            df=clusterPathways[clusterPathways$Cluster==cl,]
            if(!empty(df)){
               A = computeEdgeWeights(A,df)
               colnames(A)<-c("p1","p2","weight")
               A = aggregate(A$weight, by=list(A$p2,A$p1), FUN=sum)
            }
            #left at here: check that pathways are represented by numeric ids
            #that is strange; pathways should have had names
         }
         
      }         
      
   }

