library(data.table)
library(plyr)
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
geneSim<-function(file=bioProcessSimFile){
    message(file)
   dat<- read.csv(file,header=F,sep=",",stringsAsFactors = F)
   colnames(dat)<-c("path","genelist")
   #how many genes are there?
   genes<-unique(unlist(as.list(strsplit(dat$genelist, "\t"))))
   numGenes = length(genes)
   #how many pathways?
   numPathways = nrow(dat)
   
   vectors<-data.frame()
   for(gene in genes){
      vec<- (rownames(dat[dat$genelist %like% gene, ]))
      vectors<-rbind.fill(vectors,data.frame(gene,paste(vec,collapse=" ")))
   }
   colnames(vectors)<-c("gene","pathways")
   levels(droplevels(vectors$pathways))
   simMatrix = matrix(nrow=numGenes,ncol=numPathways)
   
   for(i in 1:numGenes){
      for(j in (i+1):numGenes){
         vec1 = strsplit(vectors[i, ]," ")
         vec2 = strsplit(vectors[j, ],"\t")
         sim = cosineSim(vec1,vec2)
         if(sim>0)
         message(i," ", j," ",sim)
      }
   }
                    
                    boolean<-false  
}



for (simFile in simFiles) {
    message(simFile)
    geneSim(simFile)
}
