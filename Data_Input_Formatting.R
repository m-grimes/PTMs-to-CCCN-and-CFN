
# Load additional packages
# library ("tidyverse")	has the following packages embedded	
#•	ggplot2, for data visualisation.
#•	dplyr, for data manipulation.
#•	tidyr, for data tidying.
#•	readr, for data import.
#•	purrr, for functional programming.
#•	tibble, for tibbles, a modern re-imagining of data frames.
#•	stringr, for strings.
#•	forcats, for factors.	

library(plyr)
library(gplots)
library(tidyverse)

#Utility Functions
"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
without <- function(x, y) x[!x %in% y] #--  x without y
nmissing <- function(x) sum(is.na(x))
filled <- function (x) {length(x) - nmissing(x)}
fractNA <- function(df) {
  result <- nmissing(df)/(dim(df)[1]*dim(df)[2])
  return(result)
}
mean.na <- function(x) mean(x, na.rm=TRUE)
max.na <- function(x) max(x, na.rm=TRUE)
min.na <- function(x) min(x, na.rm=TRUE)
sd.na <- function(x) sd(x, na.rm=TRUE)
outersect <- function(x,y){sort(c(setdiff(x,y), setdiff(y,x)))}

zero.to.NA <- function(df) {
  zer0 <- which(df==0, arr.ind = TRUE)
  cfNA <- as.matrix(df)
  cfNA[zer0] <- NA
  cfNA <- data.frame(cfNA)
  return(cfNA)
}

#########################
# Functions to import data and generate "Peptide.Name"
# requires data table with rownames 
clust.data.from.vec <- function(vec, tbl) {
  if(class(vec)=="list") {vec <- unlist(vec)}
  at <- tbl[vec,]
  acol <- names(at[,which(numcolwise(filled)(at) != 0)])
  if(length(acol)  == 1) {
    ats <- data.frame(cbind (rownames(at), as.numeric(at[, acol])))
    names(ats) <- c("Gene.Name", acol)
  }
  if(length(acol) >= 2) {
    ats <- cbind(rownames(at), at[, acol])
    names(ats)[1] <- "Gene.Name" }
  clust.data <- ats
  return (clust.data)	}
# end clust.data.from.vec

# A general version of get.gene.name; 
# use with sapply	
# Fixes some common misnamed genes and the problem that Excel turns some gene names into dates.
get.gene <- function(cell) {  
  fixgenes = c("CDC2", "2-Sep", "3-Sep", "4-Sep", "5-Sep", "7-Sep", "8-Sep", "9-Sep", "10-Sep", "11-Sep", "15-Sep", "6-Sep", "1-Oct", "2-Oct", "3-Oct", "4-Oct", "6-Oct", "7-Oct", "11-Oct", "1-Mar", "2-Mar", "3-Mar", "4-Mar", "5-Mar", "6-Mar", "7-Mar", "8-Mar", "9-Mar", "10-Mar", "11-Mar", "C11orf58", 'C17orf57', 'C3orf10',  'C7orf51', "C11orf59", "C4orf16")
  corrects = c("CDK1", "SEPT2", "SEPT3", "SEPT4", "SEPT5", "SEPT7", "SEPT8", "SEPT9", "SEPT10", "SEPT11", "SEPT15", "SEPT6", "POU2F1", "POU2F2", "POU5F1", "POU5F1", "POU3F1", "POU3F2", "POU2F3", "MARCH1", "MARCH2", "MARCH3", "MARCH4", "MARCH5", "MARCH6", "MARCH7", "MARCH8", "MARCH9", "MARCH10", "MARCH11", "SMAP", "EFCAB13", "BRK1", "NYAP1", "LAMTOR1", 'AP1AR')
  x<-unlist(strsplit(as.character(cell), ";"))	
  for (i in 1:length(fixgenes)) {
    if (x[1] == fixgenes[i]) {
      x[1] <- corrects[i] 
      return (x[1])
    } 
  }
  return(x[1]) }

# Make peptide names using this function (revised to Karen's reformatting). 
# Example:
# tencellack.head$Peptide.Name <- mapply(name.peptide, genes=tencellack.head$LeadingGeneSymbols, sites= tencellack.head$Positions, modification="ack", aa=tencellack.head$Amino.acid)
name.peptide <- function (genes, modification="p", sites, aa)	{
  genes.v <- unlist(strsplit(genes, ";", fixed = TRUE))
  sites.v <- unlist(strsplit(sites, ";", fixed = TRUE))
  sites.v <- sapply(sites.v, function (x) paste (aa, x, sep=""))
  Peptide.v <- as.character(noquote(paste(genes.v[1:length(genes.v)], modification, sites.v[1:length(sites.v)], sep=" ")))
  Peptide <- paste(unique(Peptide.v), collapse="; ")
  return(Peptide)
}

##################
# Heatmap graphing functions

hclust2 <- function(z) hclust(z, method="single")

dist2 <- function(m) {
  dm <- dist(m, method = "euclidean", diag = FALSE, upper = FALSE)
  dm[is.na(dm)] <- 2*max(abs(dm), na.rm=TRUE)
  return(dm)
}	

deblandify <- function(data.file, limit) {
  bigneg <- which(data.file<=-limit, arr.ind = TRUE)
  bigpos <- which(data.file>=limit, arr.ind = TRUE)
  data.file[bigneg] <-  -limit
  data.file[bigpos] <-  limit
  return(data.file)
}

graph.clust6d.l <- function(cluster.data.index) {
  cluster.df <- data.frame(cluster.data.index)
  if (any(grepl("Gene.Name", colnames(cluster.df)))) {
    cluster.df	<- cluster.df[, -grep("Gene.Name", colnames(cluster.df))] }
  if (ncol(cluster.df) < 2) {
    cat ("\n","This is a single sample cluster!", "\n") 
    return (cluster.df) } else {
      cluster.m <- data.matrix(cluster.df)
      #cluster.m[ is.nan (cluster.m) ] <- 0	
      #cluster.m[ is.na (cluster.m) ] <- 0
      rownames(cluster.m) <- rownames(cluster.df)
      # take out order
      cluster.mo <- cluster.m
      # cluster.mo <- cluster.m[order(-as.vector(colwise(sum.na)(data.frame(t(cluster.m))))), order(-as.vector(numcolwise(sum.na)(data.frame(cluster.m))))]
      #rbyheatcolors <- colorRampPalette(colors=c('#0000FF',  '#FFFF00'), bias=0.25, space="rgb", interpolate = "spline")
      rbyheatcolors <- colorRampPalette(colors=c('#3333FF', '#FFFF00'), bias=0.25, space="rgb", interpolate = "linear")
      # royal blue to yellow
      palette(c(rbyheatcolors(500)))
      dev.new()  # can resize manually
      heatmap.2(cluster.mo, dendrogram="row", trace="none", symbreaks=TRUE, na.color="black", labRow=rownames(cluster.mo), labCol=colnames(cluster.mo), Rowv=TRUE, Colv=NA, distfun=dist2, hclustfun=hclust, scale="none", col=palette(), colsep=NULL, rowsep=NULL, sepwidth=c(0,0), revC=FALSE, keysize=0.85, density.info='histogram', cexRow = 1, cexCol = 1, denscol="green", main=names(cluster.data.index)[1], margins=c(14, 10)) 
      return (cluster.mo) }
}

# version without dendogram and no sorting
graph.clust6nodnosort.l <- function(cluster.data.index) {
  cluster.df <- data.frame(cluster.data.index[])
  if (any(grepl("Gene.Name", colnames(cluster.df)))) {
    cluster.df <- cluster.df[, -grep("Gene.Name", colnames(cluster.df))] }
  if (ncol(cluster.df) < 2) {
    cat ("\\n","This is a single sample cluster!", "\\n") 
    return (cluster.df) } else {
      cluster.m <- data.matrix(cluster.df)
      #cluster.m[ is.nan (cluster.m) ] <- 0 
      #cluster.m[ is.na (cluster.m) ] <- 0
      rownames(cluster.m) <- rownames(cluster.df)
      #  order alphabetically
      cluster.mo <- cluster.m
      # cluster.mo <- cluster.m[order(-as.vector(colwise(sum.na)(data.frame(t(cluster.m))))), order(-as.vector(numcolwise(sum.na)(data.frame(cluster.m))))]
      #rbyheatcolors <- colorRampPalette(colors=c('#0000FF',  '#FFFF00'), bias=0.25, space="rgb", interpolate = "spline")
      rbyheatcolors <- colorRampPalette(colors=c('#3333FF', '#FFFF00'), bias=0.25, space="rgb", interpolate = "linear")
      # royal blue to yellow
      palette(c(rbyheatcolors(500)))
      dev.new()  # can resize manually
      heatmap.2(cluster.mo, dendrogram="none", trace="none", symbreaks=TRUE, na.color="black", labRow=rownames(cluster.mo), labCol=colnames(cluster.mo), Rowv=FALSE, Colv=NA, distfun=dist2, hclustfun=hclust, scale="none", col=palette(), colsep=NULL, rowsep=NULL, sepwidth=c(0,0), revC=FALSE, keysize=0.85, density.info='histogram', cexRow = 1, cexCol = 1, denscol="green", main=names(cluster.data.index)[1], margins=c(14, 10)) 
      return (cluster.mo) }
}

# Smaller font for tall graphs
graph.clust6d.la <- function(cluster.data.index) {
  cluster.df <- data.frame(cluster.data.index)
  if (any(grepl("Gene.Name", colnames(cluster.df)))) {
    cluster.df	<- cluster.df[, -grep("Gene.Name", colnames(cluster.df))] }
  if (ncol(cluster.df) < 2) {
    cat ("\n","This is a single sample cluster!", "\n") 
    return (cluster.df) } else {
      cluster.m <- data.matrix(cluster.df)
      #cluster.m[ is.nan (cluster.m) ] <- 0	
      #cluster.m[ is.na (cluster.m) ] <- 0
      rownames(cluster.m) <- rownames(cluster.df)
      # take out order
      cluster.mo <- cluster.m
      # cluster.mo <- cluster.m[order(-as.vector(colwise(sum.na)(data.frame(t(cluster.m))))), order(-as.vector(numcolwise(sum.na)(data.frame(cluster.m))))]
      #rbyheatcolors <- colorRampPalette(colors=c('#0000FF',  '#FFFF00'), bias=0.25, space="rgb", interpolate = "spline")
      rbyheatcolors <- colorRampPalette(colors=c('#3333FF', '#FFFF00'), bias=0.25, space="rgb", interpolate = "linear")
      # royal blue to yellow
      palette(c(rbyheatcolors(500)))
      dev.new()  # can resize manually
      heatmap.2(cluster.mo, dendrogram="row", trace="none", symbreaks=TRUE, na.color="black", labRow=rownames(cluster.mo), labCol=colnames(cluster.mo), Rowv=TRUE, Colv=NA, distfun=dist2, hclustfun=hclust, scale="none", col=palette(), colsep=NULL, rowsep=NULL, sepwidth=c(0,0), revC=FALSE, keysize=0.5, density.info='histogram', cexRow = 0.8, cexCol = 1, denscol="green", main=names(cluster.data.index)[1], margins=c(14, 10)) 
      return (cluster.mo) }
}
