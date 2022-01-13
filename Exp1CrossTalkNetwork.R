rm(list = ls())
library(data.table)
library(plyr)
library(tidyverse)
library(igraph)

# in RStudio: set the working directory to be directory of this R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Learn where data files are stored.
source("dataConfig.R")

# Potentially, we can compute similarity of two pathways by using different gene similarity measures
# Store the names of similarity files
simFiles  = c(bioProcessSimFile,
              molecularSimFile,
              proDomainSimFile,
              subcellularSimFile
              )

# cosineSim computes the standard cosine similarity of two gene vectors.
# Each entry in the gene vector is a pathway
# the function assumes that vec1 and vec2 do not contain duplicate values.
cosineSim <- function(vec1, vec2, len1) {
   len2 = length(vec2)
   result = 0.0
   if (is.na(vec1) || len1 == 0 ||
       is.na(vec2) || len2 == 0) {
      return(0.0)
   }
   
   
   inter = length(intersect(vec1, vec2))
   if (inter == 0) {
      result = 0.0
   }
   else{
      result = inter / (sqrt(len1) * sqrt(len2))
   }
   
   return(result)
}

# GeneSim computes similarity of two genes by analyzing the pathways they co-appear
# The simFile has pathways and genes in them.
geneSim <- function(simFile, testing = FALSE) {
   dat <- read.csv(
      simFile,
      header = FALSE,
      sep = ",",
      stringsAsFactors = F
   )
   colnames(dat) <- c("path", "genelist")
   if(length(dat$path)!=length(unique(dat$path))){
      message("Some pathways appear multiple times in ",simFile)
      return(-1);
   }
   # give ids to pathways
   numPathways <- nrow(dat)
   pathwayId <- transform(dat, id = seq_len(numPathways))
   # how many genes are there?
   genes <- sort(unique(unlist(as.list(strsplit(
      dat$genelist, "\t"
   )))))
   numGenes = length(genes)
   geneId<-transform(data.frame(genes), id = seq_len(numGenes))
   colnames(geneId) <- c("gene", "id")
   
   # how many pathways?
   if (testing) {
      message(numGenes, " genes, ", numPathways, " pathways")
   }
   
   # we will create a new 1-to-1 mapping gene->pathway
   newDat <- data.frame(gene = NA, pathway = NA)
   pathwayColumn <- numeric()
   geneColumn <- numeric()
   i <- 0
   for (row in 1:nrow(pathwayId)) {
      thisPathway <- pathwayId[row, "id"]
      for (geneName in unlist(as.list(strsplit(pathwayId[row, "genelist"], "\t")))) {
         i = i + 1
         pathwayColumn[i] <- thisPathway
         geneColumn[i] <- geneId[geneId$gene==geneName,"id"]
      }
   }
   newDat <- data.frame(geneColumn, pathwayColumn, stringsAsFactors = FALSE)
   colnames(newDat) <- c("gene", "pathway")
   # we will cache all pathways of a gene in a map to speed up future comps. 
   pathwaysOfGenes <- list()
   for (currentGene1 in 1:numGenes) {
      pathwaysOfGenes[[currentGene1]] <- newDat[newDat$gene == currentGene1,]$pathway
   }
   
   # we now start computing the pairwise similarities of genes
   geneColumn1 <- character()
   geneColumn2 <- character()
   simColumn <- numeric()
   index = 0
   message(numGenes," genes exist.")
   for(pathway1 in unique(newDat$pathway)){
      p1genes= newDat[newDat$pathway==pathway1,"gene"]
      numGenesInPathway=length(p1genes)       
      for (g1index in 1:(numGenesInPathway-1)) {
         currentGene1=p1genes[[g1index]]
         vec1 = pathwaysOfGenes[[currentGene1]]
         gene1Name=(geneId[geneId$id==currentGene1,"gene"])
         len1 <- length(vec1)
         for (g2index in (g1index + 1):numGenesInPathway) {
            
            currentGene2=p1genes[[g2index]]
            gene2Name=(geneId[geneId$id==currentGene2,"gene"])
            vec2 = pathwaysOfGenes[[currentGene2]]
            sim= cosineSim(vec1, vec2, len1)
            if(sim>0){
               index = index + 1
               geneColumn1[index] <- gene1Name
               geneColumn2[index] <- gene2Name
               simColumn[index] <- sim
            }
         }
      }
   }
   
   simDat <- data.frame(geneColumn1, geneColumn2, simColumn, stringsAsFactors = FALSE)
   return(unique(simDat))
}

loadSim<-function(simFile){
   message(simFile, " is being processed to create gene similarity values.")
   result = ""
   resFile <-  paste(
            str_sub(simFile, start = -20),
            "geneSim.rds",
            sep = ""
         )
   tryCatch({
      result <-readRDS(file = resFile)
      message("Previously computed similarity results are found and loaded.")
   }, error = function(e) {
      message("There is no previous save. Gene similarity will be computed now.")
      result = geneSim(simFile)
      saveRDS(result,
              file = resFile)
   }, finally = {
      
   })
   return(result)
}



results<-data.frame()
for(simFile in simFiles){
   geneSimMap = as.data.table(loadSim(simFile))
   colnames(geneSimMap)<-c("gene1","gene2","similarity")
   load(bioPlanetsRdata)
   
   pathways = read.csv(pathwayFile)
   for(grIndex in 1:length(pathway.net.list)){
      gr=graph_from_edgelist(as.matrix(pathway.net.list[[grIndex]][,c("source","target")]))
      # graph is here. 
      
      # How do we know that this specific crosstalk pathway network is any good?
      # We can quantify edge similarity. If two pathways are connected, 
      # we can compute their similarity by taking the average similarity of their genes.
      # 1st assumption: edges between highly similar pathways are good
      sumOfAllPathwaySim=0.0
      for(pathway1 in V(gr)$name){
         pathway1Genes <- as.character(pathways[pathways$PATHWAY_NAME==pathway1,]$GENE_SYMBOL)
         sumOfPathway1Sim<-0.0
         neList <- neighbors(gr,pathway1)$name
         df10<-geneSimMap[geneSimMap$gene1%in%pathway1Genes,]
         df20<-geneSimMap[geneSimMap$gene2%in%pathway1Genes,]
         for (pathway2 in neList){
            # strange: Bioplanet has gene symbols that are 9818 Levels: 1-Dec 1-Sep 2-Sep 4-Sep 5-Sep 6-Mar 7-Sep 9-Sep A1BG A1CF A2M ... ZYX
            pathway2Genes = as.character(pathways[pathways$PATHWAY_NAME==pathway2,]$GENE_SYMBOL)
            df1<-df10[df10$gene2%in%pathway2Genes,]
            df2<-df20[df20$gene1%in%pathway2Genes,]
            both = length(intersect(pathway1Genes,pathway2Genes))
            sumOfSim = both+sum(df1$similarity)+sum(df2$similarity)
            
            avgSimOfPathway1to2<-sumOfSim/(length(pathway1Genes)*length(pathway2Genes))
            sumOfPathway1Sim=sumOfPathway1Sim+avgSimOfPathway1to2
         }
         neighLen<-length(neList)
         if(is.finite(neighLen)&neighLen>0)
            sumOfAllPathwaySim = sumOfAllPathwaySim+sumOfPathway1Sim/neighLen
      }
      
      avgGraphSim = sumOfAllPathwaySim/length(V(gr))
      x=c(graph=grIndex,avgSim=avgGraphSim,gorder=gorder(gr),gsize=gsize(gr),simFile=simFile)
      message(grIndex," ",avgGraphSim," ",gorder=gorder(gr)," ",gsize=gsize(gr)," ",simFile)
      results= bind_rows(results,x)
   }   
}

for(sim in simFiles){
   geneSimMap = as.data.table(loadSim(sim));
colnames(geneSimMap)<-c("gene1","gene2","similarity"); 
plot2<-ggplot(data = geneSimMap,aes(x=similarity))+  
geom_histogram (binwidth = 0.05, color = "white") +
   labs (title=paste(substring(sim,first=79,99)), x="Gene similarity", y="Gene pairs")
plot2
ggsave(filename = paste(substring(sim,first=79,99),"sim.png",sep=""),plot=plot2,width=10,unit="cm")
}

results2<-results
results2$graph<-as.numeric(results2$graph)
results2$simFile<-substring(results2$simFile,first=79,99)
results2$avgSim<-round(as.numeric(results2$avgSim),3)
ggplot(data=results2,aes(x=graph,y=avgSim,group=simFile,color=simFile))+geom_line(size=2)+
   theme(legend.text=element_text(size=rel(1.3)))
results2$gorder<-as.numeric(results2$gorder)
results2$gsize<-as.numeric(results2$gsize)
ggplot(data=results2,aes(x=graph,y=gorder))+geom_line(size=1) 
ggplot(data=results2,aes(x=graph,y=gsize))+geom_line(size=1)+scale_y_log10() 
