# Protein-Protein Interaction Networks
# June 5, 2020
# Mark Grimes

# =============================================================================
# Construct a network of bioplanet pathway relationships based on lung cancer PTM data
#______________________________________________________________________________
# =======================     
# Pick up from PPI_Networks.R
#_______________________________________________________________
#
data_path <- "/Users/karenross/Documents/Bioinformatics/PIR/GrimesGroupCollab/TenPTMAnalysis/TenCellLineNetwork-Rfiles2/Rfiles/"
code_path <- "/Users/karenross/Documents/Bioinformatics/PIR/GrimesGroupCollab/PTMs-to-CCCN-and-CFN/"
source("/Users/_mark_/Dropbox/_Work/R_/MG_packages.R")
# load work
load(file=paste(data_path, "GZ_PPI_Networks2.RData", sep=""))
load(file=paste(data_path, "TenCell.RData", sep=""))
load(file=paste(data_path, "BioPlanetNetworks.RData", sep=""))
#_________________________________________#_________________________________________
#
# Calculate the number of bioplanet pathways for each gene, and the number of overlapping genes in each pathway.
# there are 1657 pathways
hist(sapply(bioplanet, length), breaks=100, col="magenta")
mean(sapply(bioplanet, length))
# [1] 44.74834
# How many genes are represented?
length(unique(unlist(bioplanet))) # 9818
length((unlist(bioplanet))) # 74148; each gene is in an average of 7.5 pathways
length(grep("GUK1", bioplanet)) # 7; 
length(grep("PAG1", bioplanet)) # 12
bioplanet[grep("PAG1", bioplanet)]
length(grep("FYN", bioplanet, fixed=T)) # 82 note that these will get partial matches
intersect(names(bioplanet[grep("PAG1", bioplanet)]), names(bioplanet[grep("FYN", bioplanet)]))
# Interesting that PDGFB signaling pathway is the only RTK pathway; others are disease, immune system and "CXCR4 signaling pathway"
# How many pathways for each gene?
# Try stringr functions
# See _work_Jun4.R for tests

bioplanetgenes <- data.frame(Gene.Name = sort(unique(unlist(bioplanet))))

length(which(str_detect("CTTN", bioplanetgenes$Gene.Name)))
str_count("FYN", test.2[[1]]) # vector
# But partial matches still a problem 
# okay try this
calculate.pathways.per.gene <- function(gene, pathwaylist) {
  sum.paths <- sum(unlist(sapply(pathwaylist, function (x) which(gene %in% x))))
  return(sum.paths)
}


# Works!
# Deal with overlapping names

# OKay go for it!
bioplanetgenes$no.pathways <- sapply(bioplanetgenes$Gene.Name, calculate.pathways.per.gene, pathwaylist=bioplanet)
# This still took over an hour!
hist(bioplanetgenes$no.pathways, breaks=100, col="aquamarine", ylim=c(0,200))
head(bioplanetgenes[order(bioplanetgenes$number.of.pathways, decreasing=TRUE),], 30)
# Behold the usual suspects. 
# Fix
bioplanetgenes <- bioplanetgenes[,c(1,3)]
#
list.common <- function (list1, list2, keeplength=3) {
  parse <- lapply(list1, function (y) sapply(list2,  function(x) intersect(x, y)))
  dims <- lapply(parse, function (x) sapply(x, length))
  keep <- which(sapply(dims, sum) > keeplength)
  pare <- parse[keep]
  prune <- lapply(pare, function (y) return (y[which(sapply(y, function (x) which(length(x) > keeplength )) > 0)]))
  newlist <- unlist(prune, recursive=FALSE)
  return(newlist)
}
bioplanet.list.common <- list.common(bioplanet, bioplanet, keeplength = 1)
# This a list of intersecting genes between all bioplanet pathways (note: includes diagonal)
matrix.common <- function(clust.list1, clust.list2)	{
  tt.dd <- matrix(data = NA, nrow = length(clust.list1), ncol = length(clust.list2), byrow = FALSE, dimnames = NULL)
  rownames(tt.dd) <- names(clust.list1)
  colnames(tt.dd) <- names(clust.list2)
  for (i in 1:length(clust.list1)) {
    for (j in 1:length(clust.list2)) {
      a <- clust.list1[[i]]
      b <- clust.list2[[j]]
      common <- intersect (a, b)
      if (length (common) >= 1 )  {
        tt.dd[i,j] <- length (common) }
      #if (length (common) >= 3 )  {	
        #cat(i, ",", j, "\t", common, "\n") 
        #if(any(fungenes %in% common))  print ("     ^-------NOTE!") }
      }
    }
  return(tt.dd)	}
  
  
  

bioplanet.matrix.common <- matrix.common(bioplanet, bioplanet)
# Note: diagonal gives the number of genes in each pathway
# Determine the number of distinct genes so that the fraction of genes in common may be calculated
matrix.outersect <- function(clust.list1, clust.list2)	{
  tt.dd <- matrix(data = NA, nrow = length(clust.list1), ncol = length(clust.list2), byrow = FALSE, dimnames = NULL)
  rownames(tt.dd) <- names(clust.list1)
  colnames(tt.dd) <- names(clust.list2)
  for (i in 1:length(clust.list1)) {
    for (j in 1:length(clust.list2)) {
      a <- clust.list1[[i]]
      b <- clust.list2[[j]]
      common <- outersect (a, b)
      if (length (common) >= 1 )  {
        tt.dd[i,j] <- length (common) }
    }
  }
  return(tt.dd)	}
bioplanet.matrix.outersect <- matrix.outersect(bioplanet, bioplanet)
# Divide these to give pathway relationships on a scale of 0 to  about 3.5
bioplanet.adj.matrix <- bioplanet.matrix.common/bioplanet.matrix.outersect
#_________________________________________________
# Treat this as an adjacency matrix
# BioPlanet Pathway Relationship = (number of genes shared between Pathwayi and Pathwayi2) / (number of distinct genes in  Pathwayi and Pathwayi2 )
#
bpc.work <- bioplanet.matrix.common
# round 2
bpc.work <- bioplanet.adj.matrix
diag(bpc.work) <- NA
hist(unlist(bpc.work), breaks=100, col="magenta", ylim=c(0, 200))
# There are some pathways with large numbers of overlapping genes.
# Use igraph functions, which means zeros must replace NA
bpc.work[is.na(bpc.work)] <- 0
bioplanet.g <- graph_from_adjacency_matrix(as.matrix(bpc.work), mode="lower", diag=FALSE, weighted="Weight")
#  
hist(edge_attr(bioplanet.g)[[1]], breaks=100, col="deepskyblue", ylim=c(0,200))
# 
# ___make an edge list file
bioplanetedges <- data.frame(as_edgelist(bioplanet.g))
names(bioplanetedges) <- c("source", "target")
bioplanetedges$Weight <- edge_attr(bioplanet.g)[[1]]
bioplanetedges$interaction <- "genes in common" 
# round 2
bioplanetedges$interaction <- "pathway relationship" 

# dim 207029      4
bioplanetedges[which(bioplanetedges$Weight <1),] # none for round 1; many small fractions in round 2
bioplanetgeneedges <- bioplanetedges
# Repeat this with the adj matrix
bioplanetreledges <- bioplanetedges
plot(bioplanetreledges$Weight~bioplanetgeneedges$Weight, pch=19, col=alpha("red", 0.2))

#_________________________________________#_________________________________________
# Apply this approach to determine evidence for pathways based on the CCCN
# Cluster Pathway Evidence = ∑( PTMs from gene(i) pathway in cluster /  total number of Pathways for gene(i)
#
# PTM cluster list is eu.sp.sed.gzallt; 
# Make a list of genes in clusters.
pep.clist <- eu.sp.sed.gzallt
# Two considerations:
# 1. We want genes with different PTMs to count towards evidence.
# 2. Do we want ambiguous PTMs (e.g. "TUBB p Y106; TUBB4B p Y106; TUBB3 p Y106; TUBB2A p Y106; TUBB2B p Y106") to count only proportionally to their ambigous representation? 
# This function will deal with lists containing ambiguous peptides separated by semicolons
get.gene.names.from.peps <- function(pepvec, pepsep="; ") {
  genevec=NULL
  for(i in 1:length(pepvec)) {
    x <- unlist(strsplit(as.character(pepvec[i]), pepsep, fixed=TRUE))
    genes <- unique(sapply(as.character(x),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
    genevec <- c(genevec, genes)
  }
  return(genevec)
}
get.gene.names.from.peps(pep.clist[["3.3.3"]])
# Take out "unique"
get.all.gene.names.from.peps <- function(pepvec, pepsep="; ") {
  genevec=NULL
  for(i in 1:length(pepvec)) {
    x <- unlist(strsplit(as.character(pepvec[i]), pepsep, fixed=TRUE))
    genes <- (sapply(as.character(x),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
    genevec <- c(genevec, genes)
  }
  return(genevec)
}
testvec <- c(pep.clist[["3.3.3"]], "RAB8A p Y6", "RAB8A p Y7")
get.gene.names.from.peps(testvec)
get.all.gene.names.from.peps(testvec) # Okay
# * Make the gene list
gene.clist <- lapply(pep.clist, get.all.gene.names.from.peps)
# Identify those list elements with ambiguous genes
get.ambiguous.genes <- function(pepvec, pepsep="; ") {
  ambigs <- pepvec[str_detect(pepvec, pepsep)]
  genevec=NULL
  for(i in 1:length(ambigs)) {
    x <- unlist(strsplit(as.character(ambigs[i]), pepsep, fixed=TRUE))
    genes <- (sapply(as.character(x),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
    genevec <- c(genevec, genes)
  }
  return(genevec)  
}
pepvec=testvec
str_detect(pepvec, pepsep)
get.ambiguous.genes(testvec)
get.all.gene.names.from.peps(testvec) %w/o% get.ambiguous.genes(testvec)
# * That might be useful. Make a list of ambiguous genes
ambig.gene.clist <- lapply(pep.clist, get.ambiguous.genes)
hist(sapply(ambig.gene.clist, length),  breaks=50, col="orange", ylim=c(0, 200))
ambig.gene.clist.length <- sapply(ambig.gene.clist, length)
summary(ambig.gene.clist.length)
# Mean 2.3; Max 147. Note that some clusters have multiple ambigous genes, so using length may not be perfect
look <- ambig.gene.clist[ambig.gene.clist.length>1]
# Make a list without ambiguous genes
gene.clist.no.ambigs <- mapply(without, x=gene.clist, y=ambig.gene.clist)
# ***
gene.clist.bioplanet.common <- matrix.common(gene.clist, bioplanet)
hist(unlist(gene.clist.bioplanet.common), breaks=100, col="magenta", ylim=c(0, 200))
# A few clusters have 10 or more genes in common. Max=30
gene.clist.bioplanet.common.no.ambigs <- matrix.common(gene.clist.no.ambigs, bioplanet)
hist(unlist(gene.clist.bioplanet.common.no.ambigs), breaks=100, col="green", ylim=c(0, 200))
# Fewer clusters have 10 or more genes in common.
# ***
ambigs.clist.bioplanet.common <- matrix.common(ambig.gene.clist, bioplanet)
hist(unlist(gene.clist.bioplanet.common.no.ambigs), breaks=100, col="blueviolet", ylim=c(0, 200))
# Count those list elements with ambiguous genes
count.ambiguous.gene.weights <- function(pepvec, pepsep="; ") {
  #ambigs <- pepvec[str_detect(pepvec, pepsep)]
  genevec=list() #Before this said "NULL" and genvec was of type numeric, not type list,. Got an error when trying to add more than one number to an element
  for(i in 1:length(pepvec)) {
    x <- unlist(strsplit(as.character(pepvec[i]), pepsep, fixed=TRUE))
    no.genes <- length(sapply(as.character(x),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
    genevec[[i]] <-  rep(1/no.genes, no.genes)
  }
  return(unlist(genevec))
}
count.ambiguous.gene.weights(testvec)
gene.clist.weights <- lapply(pep.clist, count.ambiguous.gene.weights)
# This gene.clist.weights is the list of gene clusters' weights with single PTMs weighted as 1 and ambigous PTMs weighted porportially to the number of possible matches. 
# Use this to weight gene matches
weighted.matrix.common <- function(clist, pathwaylist=bioplanet, clist.weights)	{
  tt.dd <- matrix(data = NA, nrow = length(clist), ncol = length(pathwaylist), byrow = FALSE, dimnames = NULL)
  rownames(tt.dd) <- names(clist)
  colnames(tt.dd) <- names(pathwaylist)
  for (i in 1:length(clist)) {
    for (j in 1:length(pathwaylist)) {
      a <- clist[[i]]
      b <- pathwaylist[[j]]
      common <- intersect (a, b)
      if (length (common) >= 1 )  {
        tt.dd[i,j] <- sum(clist.weights[[i]][is.element(gene.clist[[i]], pathwaylist[[j]])]) }
    }
  }
  return(tt.dd)	}
#

testvec1 <- get.all.gene.names.from.peps(pep.clist[["1.225.171"]])
testvec2 <- get.all.gene.names.from.peps(testvec)
testvec1.weights <- count.ambiguous.gene.weights(pep.clist[["1.225.171"]])
testvec2.weights <- count.ambiguous.gene.weights(testvec)
is.element(testvec1, testvec2)
intersect(testvec1, testvec2)
testvec1[is.element(testvec1, testvec2)]
testvec2[is.element(testvec2, testvec1)]
# Apply to weights
testvec1.weights[is.element(testvec1, testvec2)]
testvec2.weights[is.element(testvec2, testvec1)]
# Multiply or sum?
testvec1.weights[is.element(testvec1, testvec2)]*testvec2.weights[is.element(testvec2, testvec1)]
testvec1.weights[is.element(testvec1, testvec2)]+testvec2.weights[is.element(testvec2, testvec1)]
# Given that the probability of two independent events occuring simoultaneously is the muliplication of both, choose mulitiplication, BUT
# The bioplex genes should all have a probability of 1, so only use the weights from gene.clist

pathwaylist=head(bioplanet)
clist=head(gene.clist)
clist.weights=head(gene.clist.weights)
# go for it
gene.clist.bioplanet.weighted <- weighted.matrix.common(gene.clist, bioplanet, gene.clist.weights)
# *** gene.clist.bioplanet.weighted uses the weights from ambigous genes 
dev.new() 
hist(unlist(gene.clist.bioplanet.common), breaks=100, col="magenta", ylim=c(0, 200))
dev.new() 
hist(unlist(gene.clist.bioplanet.weighted), breaks=100, col="purple", ylim=c(0, 200))
# *** this is good. 95.5% NA
#----------------------------------------------------------------------
# Now deal with the number of pathways each gene returns
# build on gene.clist.weights <- lapply(pep.clist, count.ambiguous.gene.weights)
calculate.gene.weights.using.pathway <- function(pepvec, pepsep="; ", genepath.df) {
  genevec <- get.all.gene.names.from.peps(pepvec)
  # Add step to unweight genes that are in many pathways. Note:some genes are not in pathways!
  genesinpathways <- genevec %in% genepath.df$Gene.Name
  geneweights <- rep(1, length(genevec))
  # which(genevec[genesinpathways] == genepath.df$Gene.Name)
  # which(genepath.df$Gene.Name == unlist(genevec[genesinpathways]) )
  # genepath.df$Gene.Name %in% genevec[genesinpathways]
  # match(genevec[genesinpathways], genepath.df$Gene.Name)
  # use match because it retains the length of the vector!
  geneweights[genesinpathways] <- 1/genepath.df[match(genevec[genesinpathways], genepath.df$Gene.Name), "number.of.pathways"]
  countweights <- count.ambiguous.gene.weights(pepvec)
  weightvec <- geneweights*countweights
  return(unlist(weightvec))  
}

genepath.df=bioplanetgenes
pepvec=testvec
count.ambiguous.gene.weights(testvec)
get.all.gene.names.from.peps(testvec)
calculate.gene.weights.using.pathway(testvec, pepsep="; ", genepath.df=bioplanetgenes)
# 
gene.clist.pathway.weights <- lapply(pep.clist, calculate.gene.weights.using.pathway, pepsep="; ", genepath.df=bioplanetgenes)
# **** gene.clist.pathway.weights is the weights of genes porportional to both ambigous genes and the number of pathways per gene
dev.new() 
hist(unlist(gene.clist.pathway.weights), breaks=100, col="yellow")
hist(unlist())
# Max= 1, okay, many smaller values
# go for it****
gene.clist.bioplanet.pathway.weighted <- weighted.matrix.common(gene.clist, bioplanet, gene.clist.pathway.weights)
# ***** gene.clist.bioplanet.pathway.weighted is the 
# Cluster Pathway Evidence = ∑( PTMs from gene(i) pathway in cluster /  total number of Pathways for gene(i)
dev.new() 
hist(unlist(gene.clist.bioplanet.pathway.weighted), breaks=100, col="yellowgreen", ylim=c(0, 200))
# Now max is about 6. ****
# 
#----------------------------------------------------------------------
# Rule 5. The “large cluster penalty” rule. The larger the cluster, the less weight is given to the pathways supported by that cluster.
pep.clist.sizes <- sapply(pep.clist, length)
# divide rows of gene.clist.bioplanet.pathway.weighted
identical(rownames(gene.clist.bioplanet.pathway.clustsize.weighted), names(pep.clist.sizes))
gene.clist.bioplanet.pathway.clustsize.weighted <- gene.clist.bioplanet.pathway.weighted
identical(rownames(gene.clist.bioplanet.pathway.clustsize.weighted), names(pep.clist.sizes))
# TRUE
for (i in 1:length(pep.clist.sizes)) {
  gene.clist.bioplanet.pathway.clustsize.weighted[i,] <- gene.clist.bioplanet.pathway.clustsize.weighted[i,]/pep.clist.sizes[i]  
}
dev.new()
hist(unlist(gene.clist.bioplanet.pathway.clustsize.weighted), breaks=100, col="seagreen", ylim=c(0, 200))
# Now max is 0.75
# ***** This is now Cluster Pathway Evidence as defined by Rules 1:3 and 5
cluster.pathway.evidence <- gene.clist.bioplanet.pathway.clustsize.weighted

#----------------------------------------------------------------------

# How does this filter pathway-pathway relationships?
# Try tidyverse functions like filter_all() and across()

cpe.tbl <- tibble(cluster.pathway.evidence)
rownames(cpe.tbl) <- rownames(cluster.pathway.evidence)
# Warning message:
# Setting row names on a tibble is deprecated. 
#
path.ev <- (colwise)(sum.na)(cpe.tbl)
cpe.tbl.f1a <- filter_all(cpe.tbl, all_vars(. > 1))
cpe.tbl.f1b <-filter(cpe.tbl, across(everything(), ~ .x > 1))
identical(cpe.tbl.f1a, cpe.tbl.f1b)
# TRUE - but these appear to be empty. The filter_all function is depreciated in favor of across(). 
#..........................................................
# I don't understand tibbles yet
cpe.df <- data.frame(cluster.pathway.evidence)
# Fix first two names
names(cpe.df)[1:2] <- c("2.LTR.circle.formation", "4.1BB.dependent.immune.response")
pathway.evidence <- (colwise)(sum.na)(cpe.df)
hist(unlist(pathway.evidence), breaks=100, col="blue4", ylim=c(0, 200))
pathway.evidence[which(pathway.evidence>5)]
cluster.evidence <- rowSums(cpe.df, na.rm = T)
hist(unlist(cluster.evidence), breaks=100, col="violet", ylim=c(0, 200))
cluster.filled <- apply(cpe.df, 1, filled)
# Drop rows  and cols containing missing values
cpe.df2 <- cpe.df[apply(cpe.df, 1, function(x) !all(is.na(x))), apply(cpe.df, 2, function(x) !all(is.na(x)))]
# now 833 1488 -- only a few pruned. This also should work.
# purrr::discard(l,function(x) all(is.na(x)))
# Prune for clusters that have at least 2 pathways
cpe.df3 <- cpe.df[cluster.filled >=2,]
# dim 829 1657
# Find evidence for two pathways interacting
test <- head(cpe.df3)
test[1,] %>% keep(function(x) !is.na(x)) %>% names -> test.names
test[1,] %>% keep(function(x) !is.na(x) & x > 0.1) 
# ****
test.list <- list()
  for (i in 1:dim(test)[1]){
    test[i,] %>% keep(function(x) !is.na(x) & x > 0.1) %>% names -> test.list[[i]]
  }
# This also works
test.list <- apply(test, 1, function(y) y %>% keep(., function(x) !is.na(x) & x > 0.1) )
# ****
# Look for more than two pathways
test.list2 <- test.list[which(length(test.list)>=2)]
# Each pathway has interaction with each other pathway in this list element
test.result <- test.list2[[1]]
test.matrix <- matrix(test.result, nrow=length(test.result), ncol=length(test.result))
colnames(test.matrix) <- names(test.result)
rownames(test.matrix) <- names(test.result)
test.g <- graph_from_adjacency_matrix(test.matrix, mode="lower", weighted=TRUE, diag=F)
test.edges <- data.frame(as_edgelist(test.g))
names(test.edges) <- c("source", "target")
test.edges$Weight <- edge_attr(test.g)[[1]]
test.edges$interaction <- "cluster evidence" 
# OR
# Call edge weight the sum of evidence for interaction
test.matrix2 <- test.matrix+test.matrix
test.g2 <- graph_from_adjacency_matrix(test.matrix2, mode="lower", weighted=TRUE, diag=F)
test.edges2 <- data.frame(as_edgelist(test.g2))
names(test.edges2) <- c("source", "target")
test.edges2$Weight <- edge_attr(test.g2)[[1]]
test.edges2$interaction <- "cluster evidence" 
#  Test combining edges
test.edges3 <- rbind(test.edges, test.edges2)
ddply(test.edges3, .(source, target), numcolwise(sum))
# *****
# This will be used to combine information from different clusters.
# Create a function to do this for all the interacting pathways
create.pathway.network <- function(pathway.evidence, ev.threshold){
  ev.list <- apply(pathway.evidence, 1, function(y) y %>% keep(., function(x) !is.na(x) & x > ev.threshold) )
  if(length(ev.list) < 1) {stop("Threshold too high!")}
  # Look for more than two pathways
  ev.list.sizes <- sapply(ev.list, length)
  ev.list2 <- ev.list[which(ev.list.sizes>=2)]
  if(length(ev.list2) < 1) {stop("Threshold too high!")}
  # Create edge file for each list element
  edgelist <- list()
  for (i in 1:length(ev.list2)){
    edge.mat <- matrix(ev.list2[[i]], nrow=length(ev.list2[[i]]), ncol=length(ev.list2[[i]]))
    colnames(edge.mat) <- names(ev.list2[[i]])
    rownames(edge.mat) <- names(ev.list2[[i]])
    edge.mat2 <- edge.mat+edge.mat
    edge.g <- graph_from_adjacency_matrix(edge.mat2, mode="lower", weighted=TRUE, diag=F)
    edgefile <- data.frame(as_edgelist(edge.g))
    names(edgefile) <- c("source", "target")
    edgefile$Weight <- edge_attr(edge.g)[[1]]
    edgefile$interaction <- "cluster evidence" 
    edgelist[[i]] <- edgefile
  }
   # Combine edges
  pathway.network.1 <- do.call("rbind", edgelist)
  pathway.network <-  ddply(pathway.network.1, .(source, target, interaction), numcolwise(sum))
  return(pathway.network)
}
create.pathway.network(cluster.pathway.evidence[1:10,], 1)
# Okay, this works, generates error if threshold is too high.
# Let's try the max evidence first
sparse.net <- create.pathway.network(cluster.pathway.evidence, 0.2)
# 22 edges
next.net <- create.pathway.network(cluster.pathway.evidence, 0.1)
length(unique(c(next.net$source, next.net$target)))
# 235 edges 146 nodes
calc.net.density <- function(edgefile) {
  no.nodes <- length(unique(c(edgefile$source, edgefile$target)))
  no.edges <- dim(edgefile)[2]
  no.possible.edges <- 0.5*no.nodes*(no.nodes-1)
  net.density <- no.edges/no.possible.edges
  return(net.density)
}
# Do this iteratively to look at network dimensions as a function of threshold
pathway.net.list <- list()
density.list <- list()
thresh.vec <- seq(0.2, 0.01, -0.01)
for (i in 1:(length(thresh.vec))){
  pathway.net.list[[i]] <- create.pathway.network(cluster.pathway.evidence, thresh.vec[i])
  density.list[[i]] <- calc.net.density(pathway.net.list[[i]] )
}
#
net.edges <- sapply(pathway.net.list, function(x) dim(x)[1])
dev.new()
plot(unlist(density.list) ~ thresh.vec)
plot(net.edges ~ thresh.vec)
# Threshold values from 0.11 to 0.16 produce a sort of local maximum. 
  
  
save(bioplanet,  bioplanet.list.common, bioplanet.matrix.common, bioplanet.matrix.outersect, bioplanet.adj.matrix, bioplanetgeneedges, bioplanetreledges, bioplanetgenes, get.all.gene.names.from.peps, ambig.gene.clist, ambigs.clist.bioplanet.common, count.ambiguous.gene.weights, count.ambiguous.genes, gene.clist.bioplanet.common.no.ambigs, gene.clist.no.ambigs, get.ambiguous.genes, matrix.common, weighted.matrix.common, matrix.outersect, gene.clist.bioplanet.weighted, gene.clist.weights, pep.clist, gene.clist, calculate.gene.weights.using.pathway, gene.clist.pathway.weights, gene.clist.bioplanet.pathway.weighted, gene.clist.bioplanet.pathway.clustsize.weighted, cluster.pathway.evidence, create.pathway.network, calc.net.density, pathway.net.list, density.list, file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "BioPlanetNetworks.RData", sep=""))


##$####################paste for reuse: 
bpc.work <- bioplanet.adj.matrix
diag(bpc.work) <- NA
hist(unlist(bpc.work), breaks=100, col="magenta", ylim=c(0, 200))
# There are some pathways with large numbers of overlapping genes.
# Use igraph functions, which means zeros must replace NA
bpc.work[is.na(bpc.work)] <- 0
bioplanet.g <- graph_from_adjacency_matrix(as.matrix(bpc.work), mode="lower", diag=FALSE, weighted="Weight")
#  
hist(edge_attr(bioplanet.g)[[1]], breaks=100, col="deepskyblue", ylim=c(0,200))
# 
# ___make an edge list file
bioplanetedges <- data.frame(as_edgelist(bioplanet.g))
names(bioplanetedges) <- c("source", "target")
bioplanetedges$Weight <- edge_attr(bioplanet.g)[[1]]
bioplanetedges$interaction <- "genes in common" 
# round 2
bioplanetedges$interaction <- "pathway relationship" 

# dim 
#_________________________________________#_________________________________________
# gzallt.physical.network =	edge file of the physical interaction network: CCCN, CFN, plus negative corr edges
# The simplest way to get the node attributes is to select the genes in the entire network above and the node attribute file.
genenames <- extract.gene.names.RCy3(alkep300)
alkep300.cf <- gz.cf[gz.cf$Gene.Name %in% genenames,]
# Note: this contains PTMs. If you want only the CFN (gene nodes) use:
alkep300.gene.cf <- alkep300.cf[which(alkep300.cf$Node.ID=="gene"),]