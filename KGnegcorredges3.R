################################################################################################
# Lung cancer plus ten cell line data
# July 8, 2019
# Mark Grimes

# Examine negative correlation between different PTMs on the same protein
# source("/Users/_mark_/Dropbox/_Work/R_/MG_packages.R")
################################################################################################
# Save previous version
pepcorredges.dual.neg.v1 <- pepcorredges.dual.neg
# Prepare CCCN and CFN for plotting in Cytoscape
# Use: Allt = gzdata combined with tencell trimmed: gzallt.cor
# In large file="/Users/mark_grimes/Software\ Images/GZTenCellMatrices3.RData")
# This won't open! Memory exhausted (limit reached)
# Try previous version
# cytoscaper
load(file="/Volumes/Terra_Byte/R_Archive_2/_LINCS/GZTenCellMatrices2.RData") 
# iMac Pro
load(file="/Users/_mark_/Archive/Terra_Byte/R_Archive_2/_LINCS/GZTenCellMatrices2.RData")
# >>>>>-------> Fix to include only PTMs from >2 experiments for CCCN construction
#if(dim(gzallt.cor)[1]==9400) {
#        gzallt.cor.1 <- gzallt.cor[rownames(gzallt.cor) %in% rownames(gzdata.allt), colnames(gzallt.cor) %in% rownames(gzdata.allt)]  }
#gzallt.cor <- gzallt.cor.1
# List elements taken care of below
gzcorpack <- gzallt.cor[grep(" p ", rownames(gzallt.cor)), grep(" ack ", colnames(gzallt.cor))]	
hist(as.numeric(gzcorpack), breaks=1000, col="blue4")  
gzcorpack <- na.omit(gzcorpack)
plot(density(gzcorpack))
# no obvious pattern mean 0.0436831
gzcccnpack <- gzallt.cccn[grep(" p ", rownames(gzallt.cccn)), grep(" ack ", colnames(gzallt.cccn))]	
hist(as.numeric(gzcccnpack), breaks=1000, col="purple4")  
# no obvious pattern; mean 0.723105	
# What about dually modified proteins?
# Look at all corr edges (modified from LINCS_18r.R
#  check all corr edges
pepcorredges.g <- graph.adjacency(as.matrix(gzallt.cor), mode="lower", diag=FALSE, weighted="Weight") 
hist(edge_attr(pepcorredges.g)[[1]], breaks=1000, col="magenta", xlim=c(-1,1))  # mean 0.4380804
pepcorredges.edges <- data.frame(as_edgelist(pepcorredges.g))	
pepcorredges.edges$Weight <- edge_attr(pepcorredges.g)[[1]]
pepcorredges.edges$edgeType <- "correlation" 
pepcorredges.edges $edgeType[pepcorredges.edges$Weight<=-0.5] <- "negative correlation"
pepcorredges.edges $edgeType[pepcorredges.edges$Weight>=0.5]  <- "positive correlation"
any(is.na(pepcorredges.edges$Weight)) # True
pepcorredges.edges <- pepcorredges.edges[!is.na(pepcorredges.edges$Weight),]
 # 23890876 from 42262394
names(pepcorredges.edges)[1:2] <- c("Peptide.1", "Peptide.2")
pepgene
pepcorredges.edges$Gene.1 <- sapply(pepcorredges.edges$Peptide.1, pepgene)
pepcorredges.edges$Gene.2 <- sapply(pepcorredges.edges$Peptide.2, pepgene)
# Dim 23890876 6
pepcorredges.neg <- pepcorredges.edges[pepcorredges.edges$Weight<=-0.5,]
#  3418968
pepcorredges.dual <- pepcorredges.edges[which(pepcorredges.edges$Gene.1==pepcorredges.edges$Gene.2),]
dim(pepcorredges.dual)	# 22634
pepcorredges.dual.neg <- pepcorredges.dual[pepcorredges.dual$Weight<=-0.5,]
dim(pepcorredges.dual.neg)  # 1846
tail(pepcorredges.dual.neg[grep("CTTN", pepcorredges.dual.neg$Gene.1),])
dualmodgenes <- unique(pepcorredges.dual$Gene.1)  # 543
dualmodgenes.neg <- unique(pepcorredges.dual[pepcorredges.dual$Weight<0, "Gene.1"]) # 471
# 
head(pepcorredges.dual[order(pepcorredges.dual$Weight, decreasing=FALSE),], 20)
head(pepcorredges.dual[order(pepcorredges.dual$Weight, decreasing=TRUE),], 20)
hist(pepcorredges.dual$Weight, breaks=100, col="red")  # mean 0.2833531
# are different modifications of the same protein differently correlated?
# PHOSPHORYLATION AND ACETYLATION
dualpack1 <- pepcorredges.dual[intersect(grep(" p ", pepcorredges.dual$Peptide.1), grep(" ack ", pepcorredges.dual$Peptide.2)), ] # 0
dualpack2 <- pepcorredges.dual[intersect(grep(" ack ", pepcorredges.dual$Peptide.1), grep(" p ", pepcorredges.dual$Peptide.2)), ] # 1590
dualpack <- dualpack2
hist(dualpack$Weight, col="gold", breaks=100)  # mean 0.2892499
dualpack.genes <- extract.gene.names(dualpack) # 217
dualpack.neg <- unique(dualpack2[dualpack2$Weight<0, ]) # 199 
dualpack.neg.id <- dualpack.neg[which(dualpack.neg$Gene.1==dualpack.neg$Gene.2),] # all of them, 84, ok
hist(dualpack.neg$Weight, col="gold", breaks=100)  # 
dualpack.vneg <- unique(dualpack.neg[abs(dualpack.neg$Weight)>=0.5, ]) # 162
# NOTE: ENO1; CTTN; HSP90s...
hist(dualpack.vneg$Weight, col="gold", breaks=100)  # 
dualpack.neg.genes <- unique(dualpack2[dualpack2$Weight<0, "Gene.1"]) # 84
dualpack.vneg.genes <- unique(dualpack.vneg[, "Gene.1"]) # 77
dualpack.neg.genes[dualpack.neg.genes %in% endosome.ps]
dualpack.vvneg <- unique(dualpack.neg[abs(dualpack.neg$Weight)>0.543, ]) # 127
dualpack.vvneg.genes <- unique(dualpack.vvneg[, "Gene.1"]) # 68
# ****

dualfocus1 <- dualpack.vneg.genes[dualpack.vneg.genes %in% cytoskeletonPs]
dualfocus2 <- dualpack.vneg.genes[dualpack.vneg.genes %in% endosome.ps]
# same as 1
dualfocus <- unique(c(dualfocus1, dualfocus2))
#. 29! ****

# *****
# Graph this with neg corr edges
dualpack.vneg.edges <- dualpack.vneg[, 1:4] 
dualpack.vneg.edges$Alt.Weight <- exp (dualpack.vneg.edges$Weight*20)
names(dualpack.vneg.edges)[1:2] <- c("Gene.1", "Gene.2")
dualpackgenenames <- extract.gene.names(dualpack.vneg.edges) # 77 genes
# Need to update graphNetworkPath3
# dual.result <- graphNetworkPath3(nodenames=c("Dually", " modified proteins"), dualfocus.edges, ptmedgefile=rbind(essl.optnet.edges, dualpack.vneg.edges), datacolumns=gegzexpts, geneatts=essl.netatts, ptmcccnatts=essl.ptmcccnatts)
# win1 <- dual.result[[3]]
# win1.pos <- RCy3::getNodePosition(win1, getAllNodes(win1))
#  *** Nice figure!
# This represents dually modified proteins, their interaction with bromodomain proteins, CFN PPIs, CCCN between mods, and negative correlations between p and ack modifications <(-0.5).

# PHOSPHORYLATION AND UBIQUITINATION
dualpubi1 <- pepcorredges.dual[intersect(grep(" p ", pepcorredges.dual$Peptide.1), grep(" ubi ", pepcorredges.dual$Peptide.2)), ] # 3239
dualpubi2 <- pepcorredges.dual[intersect(grep(" ubi ", pepcorredges.dual$Peptide.1), grep(" p ", pepcorredges.dual$Peptide.2)), ] # 0
dualpubi <- dualpubi1
hist(dualpubi$Weight, col="gold", breaks=100)  # mean 0.1187465
dualpubi.genes <- extract.gene.names(dualpubi) # 109
dualpubi.neg <- unique(dualpubi[dualpubi$Weight<0, ]) # 244 edges
dualpubi.neg.id <- dualpubi.neg[which(dualpubi.neg$Gene.1==dualpubi.neg$Gene.2),] # all of them, 244, ok
hist(dualpubi.neg$Weight, col="gold", breaks=100)  # 
dualpubi.vneg <- unique(dualpubi.neg[abs(dualpubi.neg$Weight)>=0.5, ]) # 385

hist(dualpubi.vneg$Weight, col="gold", breaks=100)  # 
dualpubi.neg.genes <- unique(dualpubi[dualpubi$Weight<0, "Gene.1"]) # 161
dualpubi.vneg.genes <- unique(dualpubi.vneg[, "Gene.1"]) # 147
# NOTE: ENO1; CTTN; CTNND1, DYNLL1, HSP90AA1, YWHAZ...
dualpubi.neg.genes[dualpubi.neg.genes %in% endosome.ps] # 19
dualpubi.vvneg <- unique(dualpubi.neg[abs(dualpubi.neg$Weight)>0.543, ]) # 294
dualpubi.vvneg.genes <- unique(dualpubi.vvneg[, "Gene.1"]) # 124
# ****

# Acetylation AND UBIQUITINATION
dualackubi1 <- pepcorredges.dual[intersect(grep(" ack ", pepcorredges.dual$Peptide.1), grep(" ubi ", pepcorredges.dual$Peptide.2)), ] # 4110
dualackubi2 <- pepcorredges.dual[intersect(grep(" ubi ", pepcorredges.dual$Peptide.1), grep(" p ", pepcorredges.dual$Peptide.2)), ] # 0
dualackubi <- dualackubi1
hist(dualackubi$Weight, col="gold", breaks=100)  # mean 0.1742714
dualackubi.genes <- extract.gene.names(dualackubi) # 109
dualackubi.neg <- unique(dualackubi[dualackubi$Weight<0, ]) # 279 edges
dualackubi.neg.id <- dualackubi.neg[which(dualackubi.neg$Gene.1==dualackubi.neg$Gene.2),] # all of them, 279, ok
hist(dualackubi.neg$Weight, col="gold", breaks=100)  # 
dualackubi.vneg <- unique(dualackubi.neg[abs(dualackubi.neg$Weight)>=0.5, ]) # 577

hist(dualackubi.vneg$Weight, col="gold", breaks=100)  # 
dualackubi.neg.genes <- unique(dualackubi[dualackubi$Weight<0, "Gene.1"]) # 173
dualackubi.vneg.genes <- unique(dualackubi.vneg[, "Gene.1"]) # 164
# NOTE:  CTTN; CFL1, RAN, HSP90AA1, YWHAZ, XRCC5...
dualackubi.neg.genes[dualackubi.neg.genes %in% endosome.ps] # 9
dualackubi.vvneg <- unique(dualackubi.neg[abs(dualackubi.neg$Weight)>0.543, ]) # 441
dualackubi.vvneg.genes <- unique(dualackubi.vvneg[, "Gene.1"]) # 144
# ****
#-------------------------------------------------------------------------------------------
# COMBINE peptide edge files and gene names of interest
#  make an edge list file
gzallt.cccn.edges <- data.frame(as_edgelist(gzallt.cccn.g))  # 74001 edges
names(gzallt.cccn.edges) <- c("Peptide.1", "Peptide.2")
gzallt.cccn.edges$Weight <- edge_attr(gzallt.cccn.g)[[1]]
gzallt.cccn.edges$edgeType <- "correlation" 
gzallt.cccn.edges $edgeType[gzallt.cccn.edges$Weight<=-0.5] <- "negative correlation"
gzallt.cccn.edges $edgeType[gzallt.cccn.edges$Weight>=0.5] <- "positive correlation"
#  74001 edges 
dualmodgenes.vneg <- unique(c(dualackubi.vneg.genes, dualpubi.vneg.genes, dualpack.vneg.genes))
gzallt.cccn.edges.plus <- rbind(gzallt.cccn.edges, dualpack.vneg[,1:4], dualpubi.vneg[,1:4], dualackubi.vneg[,1:4])
genepep.edges.3 <- function(nodelist, pepkey=ld.key) {
     nodelist <- unique(nodelist)
     gpedges <- pepkey[pepkey$Gene.Name %in% nodelist, 1:2]	
     names(gpedges)[1:2] <- c("Gene.1", "Gene.2")
     gpedges$edgeType <- "peptide"
     gpedges$Weight <- 1
     gpedges$Alt.Weight <- 100
     gpedges$Directed <- FALSE
     return(unique(gpedges))
}
#-----------------------------------------------------------------------------------
#____________________________________________________________________
# Next: graph the results in Cytoscape
# Nodes	
# # Make cytoscape file
# Incorporate No.Modifications and No.Samples
# First make gene data file
gzallt.key <- data.frame(Peptide.Name= rownames(gzdata.allt), Gene.Name="")
gzallt.key$Gene.Name <- sapply(as.character(gzallt.key$Peptide.Name),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])
gzallt.key.list <- dlply(gzallt.key, "Gene.Name")
gzalltsizes <- ldply(gzallt.key.list, function(x) dim(x)[1])
names(gzalltsizes) <- c("Gene.Name", "No.Modifications") 
# focus on data from peptides in two or more experiments:
gzallt.filled <- data.frame(Peptide.Name=rownames(gzdata.allt))
gzallt.filled$No.Samples <- apply(gzdata.allt[, -grep("atio",names(gzdata.allt))], 1, filled)
gzallt.pruned <- gzallt.filled[gzallt.filled$No.Samples >= 2, ]  # 9215/9215
identical(gzallt.key$Peptide.Name, gzallt.filled$Peptide.Name)  # True okay			
gzallt.key$No.Samples <- gzallt.filled$No.Samples
gzallt.gene.key <- ddply(gzallt.key, .(Gene.Name), function(x) max(x$No.Samples))  #****
identical(gzallt.gene.key$Gene.Name, gzalltsizes$Gene.Name) # True okay	
gzallt.gene.key$No.Modifications <- gzalltsizes$No.Modifications
names(gzallt.gene.key)[2] <- "No.Samples"
# Gene data files: from gzdata.allt
identical(rownames(gzdata.allt), gzallt.key$Peptide.Name) # True okay	
gzdata.allt.df <- cbind(Gene.Name=gzallt.key$Gene.Name, gzdata.allt[, -grep("atio",names(gzdata.allt))])
# Here we want a "total signal" for each gene. The data are in log2 form so unlog to add
# This is the same as log2(sum.na(2^x))
gzalltgene.data <- ddply(gzdata.allt.df, .(Gene.Name), numcolwise(function (x) log2(sum.na(2^x))) )
rownames(gzalltgene.data) <- gzalltgene.data$Gene.Name
gzalltgene.data <- gzalltgene.data[,2:45]
# Make average ratio gene data for plotting networks
#gzallt.ave.ratios <- data.frame(
#        H3122.C.ratio = rowMeans(gzalltratios.lim[, 1:3], na.rm=TRUE),
#        H3122.PR.ratio = rowMeans(gzalltratios.lim[, 4:6], na.rm=TRUE),
#        PC9.E.ratio = rowMeans(gzalltratios.lim[, 7:9], na.rm=TRUE),
#        PC9.PR.ratio = rowMeans(gzalltratios.lim[, 10:12], na.rm=TRUE)
#)
gzallt.ave.ratios <- gzdata.allt[, grep("atio",names(gzdata.allt))]
identical(rownames(gzallt.ave.ratios), gzallt.key$Peptide.Name) # True okay	
gzallt.ave.ratios.df <- cbind(Gene.Name=gzallt.key$Gene.Name, gzallt.ave.ratios)
#****>>>>
gzalltgene.ave.ratios <- ddply(gzallt.ave.ratios.df, .(Gene.Name), function (x) numcolwise(log2(mean.na(2^x))) )
# Error
# Work around unknown error by doing this in steps
gzalltgene.ave.ratios.2 <- numcolwise(function(x) 2^x)(gzalltgene.ave.ratios)
gzalltgene.ave.ratios.2$Gene.Name <- gzalltgene.ave.ratios$Gene.Name
gzalltgene.ave.ratios.3 <- ddply(gzalltgene.ave.ratios.2, .(Gene.Name), numcolwise(mean.na) )
gzalltgene.ave.ratios.4 <- numcolwise(log2)(gzalltgene.ave.ratios.3)

rownames(gzalltgene.ave.ratios.4) <- gzalltgene.ave.ratios.3$Gene.Name
gzalltgene.ave.ratios<- gzalltgene.ave.ratios.4
#-----------------------------------------------------------------------------------
gzalltgenepep.edges <- genepep.edges.3(nodelist=gzallt.gene.key$Gene.Name, pepkey=gzallt.key)
##
# Merge edge files with common names for new RCy3
class(gzalltgene.cfn$Weight) # [1] "numeric" if not:
# gzalltgene.cfn$Weight <- as.numeric(as.character(gzalltgene.cfn$Weight))
gz.cfn <- gzalltgene.cfn
# Conform to new names in RCy3
names(gz.cfn) <- c("source", "target", "Weight", "interaction")
# gzallt.cccn.edges.plus$Weight <- as.numeric(as.character(gzallt.cccn.edges.plus$Weight))
gzallt.cccnplus <- gzallt.cccn.edges.plus[,1:4]
names(gzallt.cccnplus) <- c("source", "target", "Weight", "interaction")
gzallt.gpe <- gzalltgenepep.edges[,c(1,2,4,3)]
names(gzallt.gpe) <- c("source", "target", "Weight", "interaction")
gzallt.gpe$Weight <- as.numeric(as.character(gz.gpe$Weight))
# Unpruned
gzallt.network.all <- rbind(gz.cfn, gzallt.gpe, gzallt.cccnplus)
dim(gzallt.network.all)  # 97931     4
# Prune gpe nodes not in cccn
cccn.nodes <- unique(c(gzallt.cccnplus$source, gzallt.cccnplus$target)) # 7715
cfn.nodes <- unique(c(gz.cfn$source, gz.cfn$target))
gpe.nodes <- unique(gz.gpe$source) # 1986
gzallt.gpe.pruned <- gzallt.gpe[gzallt.gpe$source %in% cccn.nodes,]
gzallt.gpe.pruned <- gzallt.gpe.pruned[gzallt.gpe.pruned$target %in% cfn.nodes,]
gzallt.network <- rbind(gz.cfn, gzallt.gpe.pruned, gzallt.cccnplus)
dim(gzallt.network)  # [1] 95160     4
## __________________________________________
# Incorportate physcial CFNs from PPI_Networks.R
# Edgefile = gzalltgene.physical.cfn.merged
gzallt.physical.network <- rbind(gzalltgene.physical.cfn.merged, gzallt.gpe.pruned, gzallt.cccnplus)
# 86421 edges
# Saved in PPI_Networks.R as "GZ_PPI_Networks2.RData"
# __________________________________________

# # Plot whole networks!
# Function to make a file for nodes
make.anynetcf <- function(edge.df, data.file, geneatts, ptmcccnatts, func.key=func.key, use=c("total", "mean", "median", "max")) {
     cf.nodes <- unique(c(as.character(edge.df[,1]), as.character(edge.df[,2])))
     cf.data <- data.frame(data.file[rownames(data.file) %in% cf.nodes,])
     if (any(grepl("Total", names(cf.data)))) {
          cf.data <- cf.data[, -grep("Total",names(cf.data))] }	
     cf <- data.frame(Peptide.Name=rownames(cf.data))
     cf$Peptide.Name <- as.character(cf$Peptide.Name)
     cf$Gene.Name <- sapply(cf$Peptide.Name,  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])
     if(identicadim(gzallt.cccnplusl(cf$Peptide.Name, cf$Gene.Name)) {genenodes=TRUE} else {genenodes=FALSE}
     # node classification: 
     cf.functions <- func.key[func.key$Gene.Name %in% cf$Gene.Name,]
     cf <- plyr::join(cf, cf.functions, by="Gene.Name", type="left")
     # Fix any genes not in func.key
     if(any(is.na(cf))) {cf[is.na(cf)] <- "undefined"}
     # cfheadcols <- dim(cf)[2]
     # quantitiative data:
     data.class <- sapply (cf.data, class)
     #	
     if (any(use=="total")) {
          cf.Total <- rowSums(cf.data[,data.class=="numeric"], na.rm=TRUE)
          cf.Total[is.na(cf.Total)] <- 0 
          cf$Total <- as.numeric(cf.Total)	}
     if (any(use=="mean")) {
          cf.Mean <- rowMeans(cf.data[,data.class=="numeric"], na.rm=TRUE)
          cf.Mean[is.na(cf.Mean)] <- 0 
          cf$Mean <- as.numeric(cf.Mean)	}
     if (any(use=="median"))	{
          cf.Median	<- apply(cf.data[,data.class=="numeric"], 1, median, na.rm=TRUE)
          cf.Median[is.na(cf.Median)] <- 0 
          cf$Median <- as.numeric(cf.Median)	}
     if (any(use=="max"))	{
          cf.Max	<- apply(cf.data[,data.class=="numeric"], 1, function(x) {
               if (all(is.na(x))) return (NA) else {return (unique(as.numeric((x[which(abs(x)==max.na(abs(x)))]))))}				})
          if(class(cf.Max)=="list") {
               cf.Mean <- rowMeans(cf.data[,data.class=="numeric"], na.rm=TRUE)
               dupmeans <- cf.Mean[which(sapply(cf.Max, length)==2)]
               dupmax <- cf.Max[which(sapply(cf.Max, length)==2)]
               for (i in 1:length(dupmax))		{
                    if(dupmeans[[i]]>=0) dupmax[[i]] <- max.na(dupmax[[i]])
                    if(dupmeans[[i]]<0) dupmax[[i]] <- min.na(dupmax[[i]])
               } 
               cf.Max[which(sapply(cf.Max, length)==2)] <- dupmax
               cf.Max <- unlist(cf.Max)
          }	
          cf.Max[is.na(cf.Max)] <- 0 
          cf$Max <- as.numeric(cf.Max)	}
     if (genenodes==TRUE)  {
          cf=cf[,names(cf) %w/o% "Peptide.Name"]
          cf$No.Samples <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "No.Samples"]
          cf$No.Modifications <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "No.Modifications"]
          cf$ppidegree <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "ppidegree"]
          cf$ppibetween <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "ppibetween"]
          cf$norm.ppibetween <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "norm.ppibetween"]
     }
     if (genenodes==FALSE) {
          cf$Node.ID <- "gene" 
          # split between gene and peptide nodes
          cf[which(cf$Peptide.Name!=cf$Gene.Name), "Node.ID"] <- "peptide"	
          cf.list <- dlply(cf, .(Node.ID))
          cf.genes <- cf.list$gene
          cf.genes$parent <- ""
          if(length(cf.genes)>0){
               cf.genes$No.Samples <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "No.Samples"]
               cf.genes$No.Modifications <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "No.Modifications"]
               cf.genes$ppidegree <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "ppidegree"]
               cf.genes$ppibetween <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "ppibetween"]
               cf.genes$norm.ppibetween <- geneatts[geneatts $Gene.Name %in% cf$Gene.Name, "norm.ppibetween"]}
          cf.peptide <- cf.list$peptide
          cf.peptide$parent <- cf.peptide$Gene.Name
          cf.peptide$No.Samples <- ptmcccnatts[ptmcccnatts $Peptide.Name %in% cf.peptide$Peptide.Name, "No.Samples"]
          cf.peptide$No.Modifications <- 1		
          cf.peptide$ppidegree <- 0
          cf.peptide$ppibetween <- 0
          cf.peptide$pepdegree <- ptmcccnatts[ptmcccnatts $Peptide.Name %in% cf.peptide$Peptide.Name, "cccndegree"]
          cf.peptide$pepbetween <- ptmcccnatts[ptmcccnatts $Peptide.Name %in% cf.peptide$Peptide.Name, "cccnbetween"]
          # Node names need to be in first column
          # This is Peptide.Name for peptide-containing cfs
          if(length(dim(cf.genes[[1]])[1])==0) {cf.genes <- NULL}
          cf <-  rbind(cf.genes, cf.peptide) }
     if(any(is.na(cf))) {cf[is.na(cf)] <- 0}
     return(cf)
}
# Function to harmonize gene and peptide data for networks
harmonize_cfs3 <- function(pepcf, genecf) {
     genecf.new <- data.frame(Peptide.Name= genecf$Gene.Name, genecf)
     genecf.new$parent <- ""
     genecf.new$Node.ID <- "gene"
     pepcf.new <- pepcf
     pepcf.new$parent <- pepcf.new$Gene.Name
     cf <- merge(genecf.new, pepcf.new, all=TRUE)
     names(cf)[1] <- "Node"
     if(any(is.na(cf))) {cf[is.na(cf)] <- 0}
     return(cf)
}
# Use function to calculate betweenness, degree
make.netatts <- function(ig.network, ig.unfiltered, keyfile, ppinet=FALSE){
        cccndegree  <- igraph::degree(ig.network, mode="all", loops=F, normalized=F)
        cccnbetween <- igraph::betweenness(ig.network)
        net.df <- data.frame(cccndegree, cccnbetween)            
        net.df$Peptide.Name <- rownames(net.df) 
        net.df <- net.df[,c(3,1,2)]
        if(ppinet==TRUE){
                names(net.df) <- c("Gene.Name", "ppidegree", "ppibetween")
                allppibetween <- betweenness(ig.unfiltered)
                allppibetween <- allppibetween[sort(names(allppibetween))]
                allbetween.df <- data.frame(allppibetween)
                allbetween.df$Gene.Name <- rownames(allbetween.df)
                net.df <- merge(net.df, allbetween.df, all=TRUE)
                net.df$norm.ppibetween <- net.df$ppibetween/net.df$allppibetween}
        netatts.df <- merge(keyfile, net.df, all=TRUE)
        netatts.df[is.na(netatts.df)] <- 0
        return(netatts.df)
}
# CCCN: 9215 nodes
gzcccn.netatts <- make.netatts(ig.network=gzallt.cccn.g, ig.unfiltered=combined.all.ppi.g, keyfile=gzallt.key, ppinet=FALSE)
# CFN: 6188 nodes
gzgene.cfn.netatts <- make.netatts(ig.network=gzalltgene.cfn.g, ig.unfiltered=combined.all.ppi.g, keyfile=gzallt.gene.key, ppinet=TRUE)

# >>>>>>>>>>>>>>>
# CCCN node file
# Examine clusters' ratio data: tencellratios.lim.log2.trimmed plus previous data
gzdata.allt.ratios <- gzdata.allt[,names(gzdata.allt)[grep("atio",names(gzdata.allt))]]
#
gzallt.cf <- make.anynetcf(edge.df=gzallt.network, data.file=gzdata.allt[, -grep("atio",names(gzdata.allt))], geneatts=gzgene.cfn.netatts, ptmcccnatts=gzcccn.netatts, func.key=func.key, use=c("total", "median", "max"))
gzallt.allcf1 <- make.anynetcf(edge.df=gzallt.network.all, data.file=gzdata.allt[, -grep("atio",names(gzdata.allt))], geneatts=gzgene.cfn.netatts, ptmcccnatts=gzcccn.netatts, func.key=func.key, use=c("total", "median", "max"))

#________________________________________
# Include ratio data from gzallt.ave.ratios
dim(gzallt.ave.ratios) # 9215, same size as gzallt.allcf
gzallt.ave.ratios$Peptide.Name <- rownames(gzallt.ave.ratios)
# Merge
gzallt.all.cf <- merge(gzallt.allcf1, gzallt.ave.ratios, by="Peptide.Name", all=TRUE)
# cccn.nodes <- unique(c(gzallt.cccnplus$source, gzallt.cccnplus$target)) # 7715/9215
gzallt.all.cf.pruned <- gzallt.all.cf[gzallt.all.cf$Peptide.Name %in% cccn.nodes,]
# CFN node file
# Gene ratios: gzalltgene.ave.ratios
gzalltgene.cf <- make.anynetcf(edge.df=gzalltgene.cfn, data.file=gzalltgene.ave.ratios, geneatts=gzgene.cfn.netatts, ptmcccnatts=gzcccn.netatts, func.key=func.key, use=c("total", "median", "max"))
gzalltgene.ave.ratios$Gene.Name <- rownames(gzalltgene.ave.ratios)
# Merge
gzalltgene.all.cf <- merge(gzalltgene.cf, gzalltgene.ave.ratios, by="Gene.Name", all=TRUE)
cfn.nodes <- unique(c(gz.cfn$source, gz.cfn$target))  # 1986/3246
gzalltgene.all.cf.pruned <- gzalltgene.all.cf[gzalltgene.all.cf$Gene.Name %in% cfn.nodes,]
#_______________

gz.cf <- harmonize_cfs3(gzallt.all.cf, gzalltgene.all.cf)
names(gz.cf)[1] <- "id"
gz.cf.pruned <- harmonize_cfs3(gzallt.all.cf.pruned, gzalltgene.all.cf.pruned)
names(gz.cf.pruned)[1] <- "id"
# Check nodes
gznet.nodes <- unique(c(gzallt.network$source, gzallt.network$target)) # 9701
outersect(gznet.nodes, gz.cf.pruned$id) # gene nodes now 0
outersect(gznet.nodes, gz.cf$id) # many
#
gznetwork.suid <- createNetworkFromDataFrames(gz.cf, gzallt.network, title="CFN plus CCCN, All Data", collection = "Interactions")
setLayoutProperties("genemania-force-directed", list(numIterations=100, defaultSpringCoefficient=0.1, defaultSpringLength=50, minNodeMass=0.001, maxNodeMass=2000, midpointEdges=250, curveSteepness=7.0e-03, isDeterministic=1, singlePartition=0, ignoreHiddenElements=1))
layoutNetwork("genemania-force-directed")
dim(gz.cf)  #  9701   46
dim(gzallt.network) # 95160     4
dim(gzallt.physical.network) # 86421     4
# Do with physical network in which edges have been merged
gznetwork.suid <- createNetworkFromDataFrames(gz.cf, gzallt.physical.network, title="CFN plus CCCN, Physical Interactions", collection = "Interactions")
setLayoutProperties("genemania-force-directed", list(numIterations=100, defaultSpringCoefficient=0.1, defaultSpringLength=50, minNodeMass=0.001, maxNodeMass=2000, midpointEdges=250, curveSteepness=7.0e-03, isDeterministic=1, singlePartition=0, ignoreHiddenElements=1))
layoutNetwork("force-directed")
###
#___________________________________________________________________________
# Groups
# https://cytoscape-working-copy.readthedocs.io/en/latest/Creating_Networks.html#grouping-nodes
# From Alex Pico
# Alt function to just create the list of nodenames
collect.nodenames <- function() {
        nodedata <- getTableColumns("node", columns = c("id", "Gene.Name", "parent", "Node.ID"))
        nodedata[grep("gene", nodedata$Node.ID), "id"]
}

nodenames<-collect.nodenames()
# Optional smaller test case:
#nodenames<-c("CDK1",    "CRK",     "DDX5")

# Use 'group create' command directly, utlilizing column-based nodeList (see "Gene.Name:")
sapply(nodenames, function(x) {
        commandsPOST(paste0('group create groupName="',x,'"',
                            ' nodeList=Gene.Name:"',x,'"'))
})

# Note: ERRORS on 2010 iMac:
getGroupInfo("ENO2")
# NOTE: after running this double clicking collapses group; dc again expands
test=listGroups()
#RCy3::commandsPOST, HTTP Error Code: 500
#url=http://localhost:1234/v1/commands/group/list
#body={
#        "network": "SUID:931568" 
#}
#Error in commandsPOST(paste0("group list", " network=\"SUID:", net.suid,  : 
#  unrecognized (table entry): Node suid: 1130573  (table name): #SHARED_ATTRS
collapseGroup(nodenames[1])
expandGroup(nodenames[1])
#deleteGroup(nodenames[1:5])
# works on laptop; try all now
# collapseGroup(nodenames) # FAILED
# Weird error, just deleted many nodes
# Try with a loop; hung after 4, try again, then restart and try backwards
for(i in length(nodenames):1) {
       print(paste(i, nodenames[i]))
        collapseGroup(nodenames[i])
        Sys.sleep(2)
}
# Replot and try again, then set styles with the following.
# Make ratio styles for all ratios  
setNodeMapping(gz.cf)
edgeDprops.RCy32()
#setCorrEdgeAppearance(gzallt.network)     
setCorrEdgeAppearance(gzallt.physical.network)     
# setEdgeWidths.RCy32(gzallt.physical.network)  # Slow
setEdgeLineWidthMapping()
ratiocols <- names(gz.cf)[grep("atio", names(gz.cf))] %w/o% "No.Modifications"
for (i in 1:length(ratiocols)){
        plotcol <- ratiocols[i]
        style.name = paste(ratiocols[i], "Style")
        print(style.name)
        setVisualStyle("default")
        setNodeColorToRatios(plotcol)    
        copyVisualStyle('default', style.name)
        setVisualStyle(style.name)
}
#_________________________________________________________________________
# Pull threads
ep300nodes.1 <- gzallt.network[grep("EP300", gzallt.network$source),] # 266
ep300nodes.2 <- gzallt.network[grep("EP300", gzallt.network$target),] # 115
ep300.0 <- rbind(ep300nodes.1, ep300nodes.2)
ep300.cf <- gz.cf[gz.cf$id %in% unique(c(ep300nodes.1$source, ep300nodes.1$target, ep300nodes.2$source, ep300nodes.2$target)),]
gznetwork.suid <- createNetworkFromDataFrames(ep300.cf, ep300.0, title="EP300 CFN plus CCCN", collection = "Interactions")
layoutNetwork("genemania-force-directed")
# Revise functions to use source and target

extract.gene.names.RCy3 <- function (peptide.edgefile)	{
        peps <- c(peptide.edgefile[,1], peptide.edgefile[,2])
        genes <- unique(sapply(peps,  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
        return(genes) }
filter.edges.0 <- function(nodenames, edge.file) {
        if("source" %in% names(edge.file)) {
            names(edge.file)[1:2] <- c("Gene.1", "Gene.2")  }
        nodenames <-as.character(nodenames)
        a = as.character(edge.file$Gene.1)
        b = as.character(edge.file$Gene.2)
        edgefile.nodes <- unique(c(a,b))
        flub <- setdiff(edgefile.nodes, nodenames) 
        # show pruned nodes (turned off)
        # if (length(flub) >= 1) { 
        # cat("\n","\t", "The following GM names do not match ","\n","\t", flub) }	
        sel.edges <- edge.file[edge.file$Gene.1 %in% nodenames & edge.file$Gene.2%in% nodenames,]
        if("Gene.1" %in% names(sel.edges)) {
                names(sel.edges)[1:2] <- c("source", "target")  }
        if(dim(sel.edges)[1] == 0) {return(NA)} else return(sel.edges) 
}
#--------\
ep300.0.genes <- extract.gene.names(ep300.0)
ep300.cfn <- filter.edges.0(ep300.0.genes, gzallt.network)
ep300.cfn.cf <- gz.cf[gz.cf$id %in% unique(c(ep300.cfn$source, ep300.cfn$target, ep300.cfn$source, ep300.cfn$target)),]
gznetwork.suid <- createNetworkFromDataFrames(ep300.cfn.cf, ep300.cfn, title="EP300 CFN", collection = "Interactions")
layoutNetwork("genemania-force-directed")

# Refer to Figure Vignette to create functions to set node shapes, etc.
setNodeAppearance <- function(cf) {
     setBackgroundColorDefault("#949494") # grey 58
     setNodeShapeDefault( "ELLIPSE")
     setNodeColorDefault( '#F0FFFF') # azure1
     setNodeSizeDefault( 100) # for grey non-data nodes
     setNodeFontSizeDefault( 22)
     setNodeLabelColorDefault( '#000000')  # black
     setNodeBorderWidthDefault( 1.8)
     setNodeBorderColorDefault( '#888888')  # gray 
     molclasses <- c("unknown", "receptor tyrosine kinase",  "SH2 protein", "SH2-SH3 protein", "SH3 protein", "tyrosine kinase",  "SRC-family kinase",   "kinase", "phosphatase", "transcription factor", "RNA binding protein")
     #  NOTE getNodeShapes(cy) returns node shapes in random order!  Define manually 
     #	*12 for RCy2; 9 for RCy3
     # there are now 24 nodeType classes
     nodeshapes <- c("ELLIPSE","ROUND_RECTANGLE", "VEE", "VEE", "TRIANGLE", "HEXAGON", "DIAMOND", "OCTAGON", "OCTAGON", "PARALLELOGRAM", "RECTANGLE")
     setNodeSelectionColorDefault(  "#CC00FF") 
     setNodeShapeMapping ("nodeType", molclasses, nodeshapes, default.shape="ELLIPSE")
     setNodeBorderWidthMapping("nodeType", c("deacetylase","acetyltransferase","demethylase","methyltransferase","membrane protein", "receptor tyrosine kinase", "G protein-coupled receptor", "SRC-family kinase", "tyrosine kinase", "kinase", "phosphatase"), widths=c(4,12,4,12,8,16,16,12,12,12,14), 'd',default.width=4)
     if (length(cf[grep("SH2", cf$Domains), 1])>0 & !all(grep("SH2", cf$Domains) %in% which(cf$nodeType %in% molclasses))) {
          setNodeShapeBypass(cf[grep("SH2", cf$Domains) %w/o% which(cf$nodeType %in% molclasses), 1], nodeshapes[3])} 
     if (length(cf[grep("RNA", cf$nodeType), 1])>0) {
          setNodeShapeBypass(cf[grep("RNA", cf$nodeType), 1], nodeshapes[11])}
     if (length(cf[grep("transcription", cf$nodeType), 1])>0) {
          setNodeShapeBypass(cf[grep("transcription", cf$nodeType), 1], nodeshapes[10])}
     if (length(cf[grep("acetyl", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("acetyl", cf$nodeType), 1], "#FF8C00")} # darkorange
     if (length(cf[grep("methyl", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("methyl", cf$nodeType), 1], "#005CE6")} # blue
     if (length(cf[grep("membrane", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("membrane", cf$nodeType), 1], "#6600CC") # purple
          setNodeShapeBypass(cf[grep("membrane", cf$nodeType), 1], nodeshapes[2])} 
     if (length(cf[grep("kinase", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("kinase", cf$nodeType), 1], "#EE0000")} # red2
     if (length(cf[grep("phosphatase", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("phosphatase", cf$nodeType), 1], "#FFEC8B")} # lightgoldenrod1
     if (length(cf[grep("receptor", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("receptor", cf$nodeType), 1], "#BF3EFF") # darkorchid1
          setNodeShapeBypass(cf[grep("receptor", cf$nodeType), 1], nodeshapes[2])} 
     if (length(cf[grep("TM", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("TM", cf$Domains), 1], "#6600CC") # purple
          setNodeShapeBypass(cf[grep("TM", cf$Domains), 1], nodeshapes[2])} 
}
# Alex Pico writes: Try to avoid using Bypasses whenever possible. They are slow; they are not saved with a style; they stick to a particular network view. I only recommend them for a handful of things that you really must override. But you can almost always simply create a new column and make a mapping instead. 
setNodeMapping <- function(cf) {
     setBackgroundColorDefault("#949494") # grey 58
     setNodeShapeDefault("ELLIPSE")
     setNodeColorDefault("#F0FFFF") # azure1
     setNodeSizeDefault(100) # for grey non-data nodes
     setNodeFontSizeDefault( 22)
     setNodeLabelColorDefault("#000000")  # black
     setNodeBorderWidthDefault( 1.8)
     setNodeBorderColorDefault("#888888")  # gray 
     molclasses <- c("unknown", "receptor tyrosine kinase",  "SH2 protein", "SH2-SH3 protein", "SH3 protein", "tyrosine kinase",  "SRC-family kinase",   "kinase", "phosphatase", "transcription factor", "RNA binding protein")
     #  NOTE getNodeShapes(cy) returns node shapes in random order!  Define manually 
     #	*12 for RCy2; 9 for RCy3
     # there are now 24 nodeType classes
     nodeshapes <- c("ELLIPSE","ROUND_RECTANGLE", "VEE", "VEE", "TRIANGLE", "HEXAGON", "DIAMOND", "OCTAGON", "OCTAGON", "PARALLELOGRAM", "RECTANGLE")
     setNodeSelectionColorDefault(  "#CC00FF") 
     setNodeShapeMapping ("nodeType", molclasses, nodeshapes, default.shape="ELLIPSE")
     setNodeBorderWidthMapping("nodeType", c("deacetylase","acetyltransferase","demethylase","methyltransferase","membrane protein", "receptor tyrosine kinase", "G protein-coupled receptor", "SRC-family kinase", "tyrosine kinase", "kinase", "phosphatase"), widths=c(4,12,4,12,8,16,16,12,12,12,14), 'd',default.width=4)
     if (length(cf[grep("SH2", cf$Domains), 1])>0 & !all(grep("SH2", cf$Domains) %in% which(cf$nodeType %in% molclasses))) {
          setNodeShapeBypass(cf[grep("SH2", cf$Domains) %w/o% which(cf$nodeType %in% molclasses), 1], nodeshapes[3])} 
     if (length(cf[grep("RNA", cf$nodeType), 1])>0) {
          setNodeShapeBypass(cf[grep("RNA", cf$nodeType), 1], nodeshapes[11])}
     if (length(cf[grep("transcription", cf$nodeType), 1])>0) {
          setNodeShapeBypass(cf[grep("transcription", cf$nodeType), 1], nodeshapes[10])}
     if (length(cf[grep("acetyl", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("acetyl", cf$nodeType), 1], "#FF8C00")} # darkorange
     if (length(cf[grep("methyl", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("methyl", cf$nodeType), 1], "#005CE6")} # blue
     if (length(cf[grep("membrane", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("membrane", cf$nodeType), 1], "#6600CC") # purple
          setNodeShapeBypass(cf[grep("membrane", cf$nodeType), 1], nodeshapes[2])} 
     if (length(cf[grep("kinase", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("kinase", cf$nodeType), 1], "#EE0000")} # red2
     if (length(cf[grep("phosphatase", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("phosphatase", cf$nodeType), 1], "#FFEC8B")} # lightgoldenrod1
     if (length(cf[grep("receptor", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("receptor", cf$nodeType), 1], "#BF3EFF") # darkorchid1
          setNodeShapeBypass(cf[grep("receptor", cf$nodeType), 1], nodeshapes[2])} 
     if (length(cf[grep("TM", cf$nodeType), 1])>0) {
          setNodeBorderColorBypass(cf[grep("TM", cf$Domains), 1], "#6600CC") # purple
          setNodeShapeBypass(cf[grep("TM", cf$Domains), 1], nodeshapes[2])} 
}

# Function to set edge appearance
setCorrEdgeAppearance <- function(edgefile) {
     setEdgeLineWidthDefault (3)
     setEdgeColorDefault ( "#FFFFFF")  # white
     edgevalues <- getTableColumns('edge',c('Weight'))
     edgevalues['Weight']<-abs(edgevalues['Weight'])
     edgevalues['Weight']<-lapply(edgevalues['Weight'], function(x) x * 5)
     #setEdgeLineWidthBypass(edgevalues[['name']], edgevalues[['Weight']])
     names(edgevalues)<-c('Width')
     loadTableData(edgevalues, table = 'edge', table.key.column = 'SUID')
     setEdgeLineWidthMapping('Width', mapping.type = 'passthrough', style.name = 'default')
     setEdgeSelectionColorDefault ( "#FF69B4")  # hotpink
     edgecolors <- col2hex(c("red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2"))
     edgecolorsplus <- col2hex(c("deeppink", "red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2", "orangered2"))
     #  red; turquois; green; magenta; blue; violet; green;  bluegreen; black; gray; turquoiseblue; orange 
     edgeTypes <- c("pp", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of",  "controls-state-change-of", "Physical interactions", "BioPlex", "in-complex-with",  'experiments',  'database',   "Pathway", "Predicted", "Genetic interactions", "correlation", "negative correlation", "positive correlation",  'combined_score', "merged" , "intersect", "peptide", 'homology', "Shared protein domains") 
     # 22 edgeTypes            
     myarrows <- c ('Arrow', 'Arrow', 'Arrow', 'Arrow', "Arrow", 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None')
     setEdgeTargetArrowMapping( 'interaction', edgeTypes, myarrows, default.shape='None')  
     matchArrowColorToEdge('TRUE')
     setEdgeColorMapping( 'interaction', edgeTypes, edgecolors, 'd', default.color="#FFFFFF")     
     }   
# 

# Function to sent node size and color to match ratio data  
setNodeColorToRatios <- function(plotcol){
     cf <- getTableColumns('node')
     if(!(plotcol %in% getTableColumnNames('node'))){
          print (getTableColumnNames('node'))
          cat("\n","\n","\t", "Which attribute will set node size and color?")
          plotcol <- as.character(readLines(con = stdin(), n = 1))
     }
     limits <- range(cf[, plotcol])
     node.sizes     = c (135, 130, 108, 75, 35, 75, 108, 130, 135)
     #	RATIO is plotted
     #	Blue is negative: Yellow positive, Green in middle
     #		
     size.control.points = c (-100.0, -15.0, -5.0, 0.0, 5.0, 15.0, 100.0)
     color.control.points = c (-100.0, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, 100.0)
     if(limits[1] < min(size.control.points)) {
          size.control.points = c (limits[1], -15.0, -5.0, 0.0, 5.0, 15.0, 100.0)
          color.control.points = c (limits[1]-1, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, 100.0)
     }
     if(limits[2] > max(size.control.points)) {
          size.control.points = c (limits[1], -15.0, -5.0, 0.0, 5.0, 15.0, limits[2])
          color.control.points = c (limits[1]-1, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, limits[2]+1)
     }
     ratio.colors = c ('#0099FF', '#007FFF','#00BFFF', '#00CCFF', '#00FFFF', '#00EE00', '#FFFF7E', '#FFFF00', '#FFE600', '#FFD700', '#FFCC00')
     setNodeColorMapping (names(cf[plotcol]), color.control.points, ratio.colors, 'c')
     lockNodeDimensions('TRUE')
     setNodeSizeMapping (names(cf[plotcol]), size.control.points, node.sizes, 'c')
     setNodeSelectionColorDefault ( "#CC00FF") 
}
# Make ratio styles for all ratios     
ratiocols <- names(gz.cf)[grep("atio", names(gz.cf))] %w/o% "No.Modifications"
for (i in 1:length(ratiocols)){
        plotcol <- ratiocols[i]
        style.name = paste(ratiocols[i], "Style")
        setVisualStyle("default")
        setNodeColorToRatios(plotcol)    
        copyVisualStyle('default', style.name)
        setVisualStyle(style.name)
}


# Individual column method
plotcol <- "H2286_DasatinibRatio"
     style.name = 'H2286Dasatinib PTM Style'
     setVisualStyle("default")
     setNodeColorToRatios("H2286_DasatinibRatio")    
     copyVisualStyle('default', style.name)
     setVisualStyle(style.name) 
     
plotcol <- "HCC827_ErlotinibRatio"
     style.name = 'HCC827_Erlotinib PTM Style'
     setVisualStyle("default")
     setNodeColorToRatios("HCC827_ErlotinibRatio")    
     copyVisualStyle('default', style.name)
     setVisualStyle(style.name)
     
plotcol <- "norm.ppibetween"
     style.name = 'Normalized PPI betweenness'
     setVisualStyle("default")
     setNodeColorToRatios("norm.ppibetween")    
     copyVisualStyle('default', style.name)
     setVisualStyle(style.name)
     
plotcol <- "ppibetween"
     style.name = 'PPI betweenness'
     setVisualStyle("default")
     setNodeColorToRatios("ppibetween")    
     copyVisualStyle('default', style.name)
     setVisualStyle(style.name)
     
     
# Nodes that are proteins are painted according to mean ratios:
     
setNodeColorToRatios("H3122.C.ratio")    
     style.name = 'Crizotinib PROTEIN Style'
     copyVisualStyle('default', style.name)
     setVisualStyle(style.name)
     
setNodeColorToRatios("H3122.PR.ratio")    
     style.name = 'PR171 H3122 PROTEIN Style'
     copyVisualStyle('default', style.name)
     setVisualStyle(style.name)
     
setNodeColorToRatios("PC9.E.ratio")    
     style.name = 'Erlotinib PC9 PROTEIN Style'
     copyVisualStyle('default', style.name)
     setVisualStyle(style.name)
     
setNodeColorToRatios("PC9.PR.ratio")    
     style.name = 'PR171 PC9 PROTEIN Style'
     copyVisualStyle('default', style.name)
     setVisualStyle(style.name)
     
     #*****     


     
     
     #setEdgeLineWidthRule(net.w, edge.attribute.name="Weight", attribute.values=net.edges$Weight, line.widths=net.edges$Weight)



#tests
     openSession()
     edgevalues <- getTableColumns('edge',c('name','EdgeBetweenness'))
     edgevalues['EdgeBetweenness']<-abs(edgevalues['EdgeBetweenness'])
     edgevalues['EdgeBetweenness']<-lapply(edgevalues['EdgeBetweenness'], function(x) x / 2000)
     setEdgeLineWidthBypass(edgevalues[['name']], edgevalues[['EdgeBetweenness']])
     
     
     
graphNetworkPath4(nodenames, path.edges,ptmedgefile=essl.cccn.edges, datacolumns=lungcols,geneatts=essl.netatts, ptmcccnatts=essl.ptmcccnatts)

graphNetworkPath4 <- function(nodenames, path.edges, ptmedgefile, datacolumns=names(ld.fc), geneatts, ptmcccnatts){
     path.nodes <- unique(c(as.character(path.edges$Gene.1), as.character(path.edges$Gene.2)))
     # Get peptides from this network
     netpeps <- unlist(sapply(path.nodes, extract.peptides, ptmedgefile))  
     pepnodes.df <- data.frame(netpeps, pep.nodes=(sapply(netpeps, function(x) unlist(strsplit(x, " "))[1])))
     netpeps <- pepnodes.df[pepnodes.df$pep.nodes %in% path.nodes, 1]
     ptm.cccn <-filter.edges.0(netpeps, ptmedgefile)
     #nf <- list()
     #for (i in 1:(length(path.nodes)-1)) {
     #     nf[[i]] <- filter.edges.0(c(path.nodes[i], path.nodes[i+1]), ppinetwork)
     #}
     net.full <- mergeEdges(path.edges)
     net.full$Alt.Weight <- net.full$Weight
     net.gene.cf <- make.anynetcf(edge.df=net.full, data.file=ldgene.fc[, datacolumns], geneatts=geneatts, ptmcccnatts=ptmcccnatts, func.key=func.key, use=c("total", "mean", "median", "max"))
     net.pep.cf <- make.anynetcf(edge.df=ptm.cccn, data.file=ld.fc[, datacolumns], geneatts=geneatts, ptmcccnatts=ptmcccnatts, func.key=func.key, use=c("total", "mean", "median"))
     # combine node attribute cytoscape file
     net.cf <- harmonize_cfs3(pepcf=net.pep.cf, genecf=net.gene.cf)
     # make gene-peptide edges
     net.gpe <- data.frame(Gene.1=net.pep.cf$Gene.Name, Gene.2=net.pep.cf$Peptide.Name, edgeType="peptide", Weight=100, Alt.Weight=1)
     # net.gpe <- genepep.edges(ptm.cccn)
     # combine edge files
     net.full <- net.full[, c("Gene.1", "Gene.2", "Weight", "edgeType", "Alt.Weight")]
     net.edges <- rbind(net.gpe, ptm.cccn, net.full)
     if (dim(ptm.cccn)[1]>500) {
          print("This network is too large to display with PTMs, so only CFN edges will be shown.")
          net.cf <- net.gene.cf 
          net.edges <- net.full
     } 
     # Graph in RCy3
     if(dim(net.edges)[1]>49){
          graphIteratively(node.data=net.cf, edge.data=net.edges)
     } else {
          net.g <- cyPlot(net.cf, net.edges)
          net.w <- CytoscapeWindow(paste(paste(nodenames,  collapse=" to "), 1+length(getWindowList(cy))), net.g)
          displayGraph(net.w)
          layoutNetwork(net.w, "genemania-force-directed")
          #setEdgeLineWidthRule(net.w, edge.attribute.name="Weight", attribute.values=net.edges$Weight, line.widths=net.edges$Weight)
          edgeDprops(net.w)
          ratioprops(net.w, net.cf, "Total")  
          showGraphicsDetails(net.w, TRUE) 
          nodeDprops.new(net.w, net.cf)
          setEdgeWidths.log(net.w, factor=1.2)
          return(list(net.cf, net.edges, net.w))   }
}

