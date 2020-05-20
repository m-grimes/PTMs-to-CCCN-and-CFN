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
              subcellularSimFile)

# cosineSim computes the standard cosine similarity of two gene vectors.
# Each entry in the gene vector is a pathway
# the function assumes that vec1 and vec2 do not contain duplicate values.
cosineSim <- function(vec1, vec2) {
   result = 0.0
   
   if (is.na(vec1) || length(vec1) == 0 ||
       is.na(vec2) || length(vec2) == 0) {
      result = 0.0
      
   }
   if (length(vec1) != length(vec2)) {
      warning("In cosine similarity vector lengths must be the same!")
      result = -1.0
      
   }
   else {
      inter = length(intersect(vec1, vec2))
      if (inter == 0) {
         result = 0.0
      }
      else{
         result = inter / (sqrt(length(vec1)) * sqrt(length(vec2)))
      }
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
   idData <- transform(dat, id = as.numeric(factor(path)))
   # how many genes are there?
   genes <- unique(unlist(as.list(strsplit(
      dat$genelist, "\t"
   ))))
   numGenes = length(genes)
   # how many pathways?
   numPathways = nrow(dat)
   if (testing) {
      message(numGenes, " genes, ", numPathways, " pathways")
   }
   newDat <- data.frame(gene = NA, pathway = NA)
   ps <- numeric()
   gs <- character()
   i <- 0
   for (row in 1:nrow(idData)) {
      p <- idData[row, "id"]
      for (g in unlist(as.list(strsplit(idData[row, "genelist"], "\t")))) {
         i = i + 1
         ps[i] <- p
         gs[i] <- g
      }
   }
   newDat <- data.frame(gs, ps, stringsAsFactors = FALSE)
   simMatrix = matrix(nrow = numGenes, ncol = numGenes)
   
   vs <- list()
   for (i in 1:numGenes) {
      g1 <- genes[i]
      vs[[g1]] <- newDat[newDat$gs == g1,]$ps
   }
   
   l1 <- character()
   l2 <- character()
   l3 <- numeric()
   index = 0
   seen <- list()
   for (row in 1:nrow(idData)) {
      p <- idData[row, "id"]
      gList <-
         unlist(as.list(strsplit(idData[row, "genelist"], "\t")))
      for (i in 1:length(gList)) {
         g1 <- gList[i]
         vec1 = vs[[g1]]
         for (j in (i + 1):length(gList)) {
            g2 <- gList[j]
            vec2 = vs[[g2]]
            index = index + 1
            l1[index] <- g1
            l2[index] <- g2
            l3[index] <- cosineSim(vec1, vec2)
         }
      }
   }
   simDat <- data.frame(l1, l2, l3, stringsAsFactors = FALSE)
   return(unique(simDat))
}


addEdges <- function(B, dic, df, testing = F) {
   pathways = sort(unique(as.character(df$PATHWAY_ID)))
   len <- length(pathways)
   for (i in 1:(len - 1)) {
      for (j in (i + 1):len) {
         res = 0
         if (is.na(B[dic[pathways[i]], dic[pathways[j]]]))
            res = 0
         else
            res = B[dic[pathways[i]], dic[pathways[j]]]
         B[dic[pathways[i]], dic[pathways[j]]] = 1 + res
      }
   }
   return(B)
}


pathwaySim<-function(pathway1,pathway2,geneSimMap,pathways){
   genes1 = as.character(pathways[pathways$PATHWAY_ID==pathway1,]$GENE_SYMBOL)
   genes2 = as.character(pathways[pathways$PATHWAY_ID==pathway2,]$GENE_SYMBOL)
   
   sim = 0.0
   for(g1 in genes1){
      d1<- geneSimMap[geneSimMap$gene1==g1,]
      d2<- geneSimMap[geneSimMap$gene2==g1,]
      for(g2 in genes2){
         if(g1==g2){
            sim=sim+1.0
         }
         else{
            s1=d1[d1$gene2==g2,][1]$similarity
            if(!is.na(s1)) sim = sim+s1
            else{ s2=d2[d2$gene1==g2,][1]$similarity
            if(!is.na(s2)) sim = sim+s2
            }
         }
         
      }
   }
   return(sim/(length(genes1)*length(genes2)))
}

loadSim<-function(simFile){
   message(simFile, " is being processed to create gene similarity values.")
   result = ""
   tryCatch({
      sim1 <-
         readRDS(file = paste(
            "results/",
            str_sub(simFile, start = -20),
            "geneSim.rds",
            sep = ""
         ))
      message("Previously computed similarity results are found and loaded.")
      result = sim1
   }, error = function(e) {
      message("There is no previous save. Gene similarity will be computed now.")
      result = geneSim(simFile)
      saveRDS(result,
              file = paste(
                 "results/",
                 str_sub(simFile, start = -20),
                 "geneSim.rds",
                 sep = ""
              ))
   }, finally = {
      
   })
   return(result)
}
simFile= bioProcessSimFile
geneSimMap = as.data.table(loadSim(simFile))
# rename the similarity columns
colnames(geneSimMap) <- c("gene1", "gene2", "similarity")
# load pathways and their genes
pathways = read.csv(pathwayFile)
# how many genes are found in the pathways file?
a <- (c(unique(geneSimMap$gene1), unique(geneSimMap$gene2)))
b <- unique(pathways$GENE_SYMBOL)
c <- intersect(a, b)
message(length(c),
        " genes from ",
        pathwayFile,
        " are also found in the similarity file ",
        simFile)

# load clusters and their genes
clusters = read.csv(clusterFile, sep = "\t", stringsAsFactors = FALSE)
# some rows contain multiple sites (e.g., TUBB4B ubi K379; TUBB2A ubi K379; TUBB2B ubi K379	1.225.171)
# we will split those and add as rows
s <- strsplit(clusters$Site, split = as.character(';'))
clusters = data.frame(V1 = unlist(s), V2 = rep(clusters$Cluster, sapply(s, length)))
colnames(clusters) <- c("Site", "Cluster")
# it is very easy to hate R
clusters$Site <- trimws(as.character(clusters$Site))
# we are interested in the symbol of the Site
# for example, we need to extract "TUBB4B" from TUBB4B ubi K379
clusters$Site = sapply(strsplit(clusters$Site, " "), `[`, 1)
allPtmCounts = data.frame(table(clusters$Cluster))
colnames(allPtmCounts) = c("Cluster", "ptmCount")
# now, we need to create pathways of clusters
clusterPathways <-
   merge(clusters,
         pathways,
         by,
         by.x = "Site",
         by.y = "GENE_SYMBOL",
         sort = TRUE)
clusterPathways$Cluster <- as.character(clusterPathways$Cluster)
# how many ptms there are in each cluster
invPtmCounts = data.frame(table(clusters$Cluster))
colnames(invPtmCounts) = c("Cluster", "ptmCount")

# how many total pathways
pNames <- unique(clusterPathways$PATHWAY_ID)
pwCount = length(pNames)
dic <- 1:pwCount
names(dic) <- pNames

# how many pathways are in each cluster
invPwCounts = data.frame(table(clusterPathways$Cluster))
colnames(invPwCounts) = c("Cluster", "pathwayCount")

# finally, we are ready for computations
maxStep = 10


percOfPathways = 2

maxInvolvedPTMs = 30000

for (step in seq.int(1, maxStep, 1)) {
   percOfPtms = step / maxStep
   
   # create a frame to hold edges
   A = data.frame(
      from = character(),
      to = character(),
      weight = numeric(),
      stringsAsFactors = FALSE
   )
   B = matrix(0L, nrow = pwCount, ncol = pwCount)
   s = 0
   # We will take each cluster, find its PTMs
   for (clID in unique(clusters$Cluster)) {
      s = s + 1
      # Just for seeing the progress
      if (s %% 100 == 0) {
        # message(s,"th cluster was processed.")
      }
      # find PTMs of the cluster
      ptmCount = allPtmCounts[allPtmCounts$Cluster == clID,]$ptmCount
      
      # rule 0: we will ignore very large clusters
      if (ptmCount >= maxInvolvedPTMs) {
         next
      }
      # Involved clusters PTMs are those who are connected to a pathway
      involvedPtmCount = invPtmCounts[invPtmCounts$Cluster == clID,]$ptmCount
      
      # how many pathways are connected to this cluster?
      involvedPathwayCount = 0
      p <- invPwCounts[invPwCounts$Cluster == clID,]
      if (nrow(p) > 0) {
         involvedPathwayCount = p$pathwayCount
      }
      
      # rule 1: how much percent of the cluster ptms are found in the pathway file
      if (ptmCount == 0 ||
          involvedPtmCount / (1.0 * ptmCount) < percOfPtms) {
         next
      }
      # rule 2: how many pathways are included per involved ptm, higher the better?
      if (involvedPathwayCount == 0 ||
          (involvedPathwayCount / involvedPtmCount) < percOfPathways) {
         next
      }
      
      # So this cluster survives the rules 1 and 2
      # we will take its ptms, find their pathways and
      # add an edge between two pathways because of this cluster.
      # each cluster adds new edges between pathways.
      df = clusterPathways[clusterPathways$Cluster == clID,]
      if (!empty(df)) {
         # B is a matrix that cumulatively records how many clusters connect two pathways
         B = addEdges(B, dic, df)
      }
   }# computations for rules 1 and 2 end here.
   
   # At this point, certain pathways have edges between them.
   # Each edge indicates that two pathways are connected through PTMs in a cluster.
   # rule 3: filtering pathway pairs by the number of common clusters
   for (minPtmClusters in seq.int(10, 200, 10)) {
      # create a copy of the data frame
      B2 = B
      # if the rule is not satisfied, remove the edge between pathways
      B2[B2 < minPtmClusters] = 0
      # rule is satisfied, set the edge to 1
      B2[B2 >= minPtmClusters] = 1
      
      # does the graph have any edge?
      if (sum(B2) > 0) {
         gr = graph_from_adjacency_matrix(B2, mode = "undirected")
         # graph is here. Now what to do with this?
         gr = simplify(gr, remove.multiple = T)
         gr = delete.vertices(simplify(gr), degree(gr) == 0)
         #plot(gr, main = paste(minPtmClusters))
         
         # How do we know that this specific crosstalk pathway network is any good?
         # We can quantify edge similarity. If two pathways are connected, 
         # we can compute their similarity by taking the average similarity of their genes.
         # 1st assumption: edges between highly similar pathways are good
         avgEdgeSim=0.0
         for(node in V(gr)){
            pathway1 = (names(dic[node]))
           
            for ( neig in neighbors(gr,node)){
               pathway2 = (names(dic[neig]))
               # strange: Bioplanet has gene symbols that are 9818 Levels: 1-Dec 1-Sep 2-Sep 4-Sep 5-Sep 6-Mar 7-Sep 9-Sep A1BG A1CF A2M ... ZYX
               
               #avgEdgeSim = sim+pathwaySim(pathway1,pathway2,geneSimMap,(pathways))
            }
         }
         message(step,",",maxInvolvedPTMs,",",percOfPathways,",",minPtmClusters,",",avgEdgeSim,",",gorder(gr),",",gsize(gr))
      }
   }
}
