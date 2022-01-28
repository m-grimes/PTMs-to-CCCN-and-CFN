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
localpath = getwd()	
comp_path <- unlist(strsplit(localpath, "Dropbox"))[[1]]
if(length(grep("Dropbox", comp_path))==0) { 
  comp_path <- unlist(strsplit(localpath, "Documents"))[[1]]}
source(paste(comp_path, "/Dropbox/_Work/R_/MG_packages.R", sep=""))
# load work
load(file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "GZ_PPI_Networks2.RData", sep=""))
load(file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "TenCell.RData", sep=""))
load(file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "BioPlanetNetworks.RData", sep=""))


#_________________________________________
#_________________________________________
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
# But partial matches still a problem   
# okay try this
calculate.pathways.per.gene <- function(gene, pathwaylist) {
  sum.paths <- sum(unlist(sapply(pathwaylist, function (x) which(gene %in% x))))
  return(sum.paths)
}

# Works!
# OKay go for it!
bioplanetgenes$no.pathways <- sapply(bioplanetgenes$Gene.Name, calculate.pathways.per.gene, pathwaylist=bioplanet)
# This still took over an hour!
hist(as.numeric(bioplanetgenes$no.pathways), breaks=100, col="aquamarine")
head(bioplanetgenes[order(bioplanetgenes$number.of.pathways, decreasing=TRUE),], 30)
# Behold the usual suspects. 
# Use correct columns (after tests)
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
      if (length (common) >= 3 )  {	
        cat(i, ",", j, "\t", common, "\n") 
        if(any(fungenes %in% common)) print ("     ^-------NOTE!")	}
    }
  }
  return(tt.dd)	}
bioplanet.matrix.common <- matrix.common(bioplanet, bioplanet)
# 1657 x 1657 pathways
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
      outer.sect <- outersect (a, b)
      if (length (outer.sect) >= 1 )  {
        tt.dd[i,j] <- length (outer.sect) }
    }
  }
  return(tt.dd)	}
bioplanet.matrix.outersect <- matrix.outersect(bioplanet, bioplanet)
# Divide these to give pathway relationships on a scale of 0 to  about 3.5
bioplanet.adj.matrix <- bioplanet.matrix.common/bioplanet.matrix.outersect
# Modifiy function to use Jaccard distance
## Jacard formula: #common / (#i + #j - #common) or length(intersect(i,j)/length(union(i,j))
matrix.union <- function(clust.list1, clust.list2)	{
  tt.dd <- matrix(data = NA, nrow = length(clust.list1), ncol = length(clust.list2), byrow = FALSE, dimnames = NULL)
  rownames(tt.dd) <- names(clust.list1)
  colnames(tt.dd) <- names(clust.list2)
  for (i in 1:length(clust.list1)) {
    for (j in 1:length(clust.list2)) {
      a <- clust.list1[[i]]
      b <- clust.list2[[j]]
      denominator <- union(a, b)
      if (length (denominator) >= 1 )  {
        tt.dd[i,j] <- length (denominator) }
    }
  }
  return(tt.dd)	}
bioplanet.matrix.union <- matrix.union(bioplanet, bioplanet)
bioplanet.jaccard.matrix <- bioplanet.matrix.common/bioplanet.matrix.union
range(unlist(bioplanet.jaccard.matrix), na.rm=T)
# 0.0005350455 to 1.0000000000 (including diagnal)
#_________________________________________________
# Treat this as an adjacency matrix
# BioPlanet Pathway Relationship = (number of genes shared between Pathwayi and Pathwayi2) / (number of distinct genes in  Pathwayi and Pathwayi2 )
# round 1:
bpc.work <- bioplanet.matrix.common
# round 2:
bpc.work <- bioplanet.adj.matrix
# round 3:
bpc.work <- bioplanet.jaccard.matrix
diag(bpc.work) <- NA
hist(unlist(bpc.work), breaks=100, col="magenta", ylim=c(0, 200))
# There are some pathways with large numbers of overlapping genes.
# Use igraph functions, which means zeros must replace NA
bpc.work[is.na(bpc.work)] <- 0
bioplanet.g <- graph_from_adjacency_matrix(as.matrix(bpc.work), mode="lower", diag=FALSE, weighted="Weight")
#  
hist(edge_attr(bioplanet.g)[[1]], breaks=100, col="deepskyblue", ylim=c(0,200))
range(edge_attr(bioplanet.g)[[1]])
# [1] 0.0005350455 0.7777777778
# ___make an edge list file
bioplanetedges <- data.frame(as_edgelist(bioplanet.g))
names(bioplanetedges) <- c("source", "target")
bioplanetedges$Weight <- edge_attr(bioplanet.g)[[1]]
bioplanetedges$interaction <- "genes in common" 
# round 2
bioplanetedges$interaction <- "pathway relationship" 
# round 3
bioplanetedges$interaction <- "pathway Jaccard similarity" 
# dim 207029      4
bioplanetedges[which(bioplanetedges$Weight <1),] # none for round 1; many small fractions in round 2
bioplanetgeneedges <- bioplanetedges
# Repeat this with the adj matrix
bioplanetreledges <- bioplanetedges
# and with jaccard similarity
bioplanetjaccardedges <- bioplanetedges
# **** Max 0.7778, similar to cluster.pathway evidence below
plot(bioplanetreledges$Weight~bioplanetgeneedges$Weight, pch=19, col=alpha("red", 0.2))
plot(bioplanetreledges$Weight~bioplanetjaccardedges$Weight, pch=19, col=alpha("forestgreen", 0.2))
# *** interesting. Max 0.7778, similar to cluster.pathway evidence below
# 
#_________________________________________
#_________________________________________
# Apply this approach to determine evidence for pathways based on the CCCN
# Cluster Pathway Evidence = ∑( PTMs from gene(i) pathway in cluster /  total number of Pathways for gene(i)
#
# PTM cluster list is eu.sp.sed.gzallt; Revised renamed> essgzallt
length(essgzallt.data) #818
length(essgzallt) #818
# Make a list of genes in clusters.
pep.clist <- essgzallt.data
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
hist(unlist(gene.clist.bioplanet.common), breaks=40, col="magenta", ylim=c(0, 200))
# A few clusters have 10 or more genes in common. Max=30
gene.clist.bioplanet.common.no.ambigs <- matrix.common(gene.clist.no.ambigs, bioplanet)
hist(unlist(gene.clist.bioplanet.common.no.ambigs), breaks=100, col="green", ylim=c(0, 200))
# Fewer clusters have 10 or more genes in common.
# ***
ambigs.clist.bioplanet.common <- matrix.common(ambig.gene.clist, bioplanet)
hist(unlist(gene.clist.bioplanet.common.no.ambigs), breaks=40, col="blueviolet", ylim=c(0, 200))
# Count those list elements with ambiguous genes
count.ambiguous.gene.weights <- function(pepvec, pepsep="; ") {
  #ambigs <- pepvec[str_detect(pepvec, pepsep)]
  genevec=NULL
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
        tt.dd[i,j] <- sum(clist.weights[[i]][is.element(clist[[i]], pathwaylist[[j]])], na.rm=TRUE) }
    }}
  return(tt.dd)	}
#

testvec1 <- get.all.gene.names.from.peps(c(pep.clist[["1.225.171"]], Rabp8="RAB8A p 8"))
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
# Sum
testvec1.weights[is.element(testvec1, testvec2)]*testvec2.weights[is.element(testvec2, testvec1)]
testvec1.weights[is.element(testvec1, testvec2)]+testvec2.weights[is.element(testvec2, testvec1)]
# test function
clist=testvec2
pathwaylist=testvec1
clist.weights=testvec2.weights
weighted.matrix.common(testvec2, testvec1, testvec2.weights)

# Given that the probability of two independent events occuring simoultaneously is the muliplication of both, choose mulitiplication, BUT
# The bioplex genes should all have a probability of 1, so only use the weights from gene.clist
# # > > > . NOTE  Test intersection to return 2 for RAB8A: 
# OKAY: Note that is.element() is used to retrieve weights, not intersect(), which is used just as a tag to get past if()
pathwaylist=head(bioplanet)
clist=head(gene.clist)
clist.weights=head(gene.clist.weights)
weighted.matrix.common(clist, pathwaylist, clist.weights)
# go for it
gene.clist.bioplanet.weighted <- weighted.matrix.common(gene.clist, bioplanet, gene.clist.weights)
# *** gene.clist.bioplanet.weighted uses the weights from ambigous genes 
dev.new() 
hist(unlist(gene.clist.bioplanet.common), breaks=100, col="magenta", ylim=c(0, 200))
dev.new() 
hist(unlist(gene.clist.bioplanet.weighted), breaks=100, col="purple", ylim=c(0, 200))
# *** this is good. 
max.na(gene.clist.bioplanet.weighted)  # 28.43056
fractNA(gene.clist.bioplanet.weighted)
# 95.5% NA
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

# Visualize?
# @# Crash@! graph.clust5(cluster.pathway.evidence)


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
# What is the size of the network without any threshold?
total.net <- create.pathway.network(cluster.pathway.evidence, 0)
#  dim (total.net) 645709      4
no.nodes.total.net <- length(unique(c(total.net$source, total.net$target))) # 1488
# no. possible edges = 0.5*no.nodes*(no.nodes-1)
no.poss.total.edges <- 0.5*(no.nodes.total.net*(no.nodes.total.net-1)) # 1106328
total.net.density <- dim (total.net)[1]/no.poss.total.edges # 0.5836506
total.net <- total.net[order(total.net$Weight, decreasing=TRUE),]
total.net.fpe <- filter.pathway.edges(total.net, bioplanetjaccardedges)
total.pathway.net <- total.net.fpe # for saving. 
names(total.pathway.net) <- names(big.path.net)
total.pathway.net <- total.pathway.net[order(total.pathway.net$Weight.clust,decreasing=TRUE),]
total.pathway.net$Combined.Weight <- total.pathway.net$Weight.bp + total.pathway.net$Weight.clust
total.pathway.net.no.bp <- total.pathway.net[which(total.pathway.net$Weight.bp==0),]
# Change the name of total.pathway.net$Weight to Weight.normalized
names(total.pathway.net)[7] <- "Weight.normalized"
# >>>>***
# Create pathway crosstalk network with individual cluster and bioplanet edges
pathway.crosstalk.network <- rbind(total.net, bioplanetjaccardedges)
# 852738 edges
##############################################################################
# 
##############################################################################
# Network characterization
calc.net.density <- function(edgefile) {
  no.nodes <- length(unique(c(edgefile$source, edgefile$target)))
  no.edges <- dim(edgefile)[2]
  no.possible.edges <- 0.5*no.nodes*(no.nodes-1)
  net.density <- no.edges/no.possible.edges
  return(net.density)
}
# Do this iteratively to look at network dimensions as a function of threshold
pathway.net.list2 <- list()
density.list2 <- list()
thresh.vec <- seq(0.2, 0.01, -0.01)
thresh.vec2 <- c(seq(0.2, 0.082, -0.002), seq(0.08, 0.01, -0.01))
for (i in 1:(length(thresh.vec2))){
  pathway.net.list2[[i]] <- create.pathway.network(cluster.pathway.evidence, thresh.vec2[i])
  density.list2[[i]] <- calc.net.density(pathway.net.list2[[i]] )
}
#
net.edges <- sapply(pathway.net.list2, function(x) dim(x)[1])
net.nodes <- sapply(pathway.net.list2, function(x) length(unique(c(x[,2], x[ ,1]))))
net.atts.df <- data.frame(threshold=thresh.vec2, nodes=net.nodes, edges=net.edges, density=unlist(density.list2))
dev.new()
plot(unlist(density.list) ~ thresh.vec)
plot(net.edges ~ thresh.vec2)
plot(net.nodes ~ thresh.vec2)
plot(log2(net.nodes) ~ thresh.vec)
# Threshold values from 0.11 to 0.16 produce a sort of local maximum. 
###>>>>>----
# Calculate background from bioplanet
test.edges <- pathway.net.list[[1]]
test.g <- graph_from_data_frame(test.edges)
test.i1 <- intersect(test.edges[,1:2], bioplanetjaccardedges[,1:2])
# reverse
test.i2 <- intersect(test.edges[,1:2], bioplanetjaccardedges[,c(2,1)])
identical(test.i1, test.i2)
# [1] TRUE
identical(test.i1, test.edges[,1:2])
# FALSE - pulls in 19 from bioplanet vs 22 edges in test.edges
# **** This means that there is a set of pathway edges not predicted, which is good!
dim(unique(intersect(test.i1, test.i2))) #19
dim(bioplanetjaccardedges[intersect(test.edges[,1:2], bioplanetjaccardedges[,1:2])],)
# Doesn't work
test.i4 <- bioplanetjaccardedges %>% right_join(test.edges, by=c("source","target"))
# Works! 
# Check for reverse
test.edges.rev <- test.edges[, c(2,1,3,4)]
names(test.edges.rev) <- names(test.edges)
test.i4r <- bioplanetjaccardedges %>% right_join(test.edges.rev, by=c("source","target"))
apply(test.i4, 1, function (x) any(is.na(x)))
apply(test.i4r, 1, function (x) any(is.na(x)))
# no rows in reverse are in Jaccard 
any(!is.na(test.i4r$interaction.x))
# I wonder if this is generally true?
 # No - 
# Attempt to recover weight to normalize 
filter.pathway.edges <- function(query.edges, edge.file=bioplanetjaccardedges) {
  sel.edges.forward <- edge.file %>% right_join(query.edges, by=c("source","target"))
  # Test for reverse edges
    query.edges.rev <- query.edges[, c(2,1,3,4)]
    names(test.edges.rev) <- names(query.edges) 
    sel.edges.rev <- edge.file %>% right_join(query.edges.rev, by=c("source","target"))
    if (any(!is.na(sel.edges.rev$interaction.x))) {
      sel.edges.rev <- sel.edges.rev[!is.na(sel.edges.rev$interaction.x),]
      sel.edges <- rbind(sel.edges.forward, sel.edges.rev)
    } else {sel.edges = sel.edges.forward}
  if(dim(sel.edges)[1] == 0) {return(NA)} else {
    # convert NA to zero
    sel.edges[is.na(sel.edges)] <- 0
    # Calculate normalized interaction
    sel.edges$Weight <- sel.edges$Weight.y - sel.edges$Weight.x
    sel.edges$interaction <- "normalized cluster evidence"
    return(unique(sel.edges)) }
}
test.i4 <- filter.pathway.edges(test.edges, bioplanetjaccardedges)
# THis works! ***
# Apply to list of networks.
pathway.net.list3 <- lapply(pathway.net.list2, filter.pathway.edges)
# Calculate attributes to check
net.edges3 <- sapply(pathway.net.list3, function(x) dim(x)[1])
net.nodes3 <- sapply(pathway.net.list3, function(x) length(unique(c(x[,2], x[ ,1]))))
net.atts3.df <- data.frame(threshold=thresh.vec2, nodes=net.nodes3, edges=net.edges3, density=unlist(density.list2))
identical(net.atts.df, net.atts3.df)  # TRUE
# Use the biggest network to examine the relationships between normalized Weight
big.path.net <- pathway.net.list3[[68]]
#plot(big.path.net$Weight.x ~ big.path.net$Weight.y)
#summary(lm(big.path.net$Weight.x ~ big.path.net$Weight.y))
# Rename
names(big.path.net) <- c("source", "target" , "Weight.bp", "bp.interaction", "clust.interaction", "Weight.clust", "Weight", "interaction" )
#
plot(big.path.net$Weight.bp ~ big.path.net$Weight.clust)
summary(lm(big.path.net$Weight.bp ~ big.path.net$Weight.clust))# R-squared:  0.01322 
# **** this is good, means that we are NOT just recapitulating relationships among pathways
plot(big.path.net$Weight ~ big.path.net$Weight.clust)
summary(lm(big.path.net$Weight ~ big.path.net$Weight.clust))
# R-squared:  0.9693; the normalization didn't move the points very much.


big.path.net <- big.path.net[order(big.path.net$Weight, decreasing=TRUE),]
hist(big.path.net$Weight, breaks=100, col="turquoise", ylim=c(0,200))
dim(big.path.net[big.path.net$Weight>4,])
# 90
new.path.net <- big.path.net[big.path.net$Weight.bp==0,]
hist(new.path.net$Weight, breaks=100, col="purple", ylim=c(0,200))
dim(new.path.net[new.path.net$Weight>1.5,])
# look at bottom end
plot(big.path.net$Weight ~ big.path.net$Weight.clust, xlim=c(-1,2))
plot(big.path.net$Weight.clust ~ big.path.net$Weight, xlim=c(-0.7,2), ylim=c(-0.1,2), pch=19, col=alpha("purple", 0.2))
# Look at smaller network 46 with Threshold 0.110: 140 nodes  217 edges
small.path.net <- pathway.net.list3[[46]]
small.path.net <- small.path.net[order(small.path.net$Weight, decreasing=TRUE),]
names(small.path.net) <- c("source", "target" , "Weight.bp", "bp.interaction", "clust.interaction", "Weight.clust", "Weight", "interaction" )
plot(small.path.net$Weight.clust ~ small.path.net$Weight, xlim=c(-0.7,2), ylim=c(-0.1,2), pch=19, col=alpha("purple", 0.2))
plot(small.path.net$Weight.clust ~ small.path.net$Weight.bp, pch=16, col=alpha("magenta3", 0.5))
summary(lm(small.path.net$Weight.clust ~ small.path.net$Weight.bp))
#  R-squared:  -0.004617
smaller.path.net <- small.path.net[small.path.net$Weight>0.2, ] # now 157 edges
plot(smaller.path.net$Weight.clust ~ smaller.path.net$Weight, pch=16, col=alpha("darkgreen", 0.5))
# Look at the top pathway-pathway interactions in networks of different thresholds

pathway.net.list4 <- lapply(pathway.net.list3, function(x) x[order(x$Weight, decreasing = TRUE),])
newnames <- c("source", "target" , "Weight.bp", "bp.interaction", "clust.interaction", "Weight.clust", "Weight", "interaction" )
pathway.net.list4 <- lapply(pathway.net.list4, setNames, newnames)
# Look at the top interactions
lapply(pathway.net.list4, head)
# interesting, the weights are different, and the top interactions change 
# due to inclusion of more evidence when chainging the threshold for individual interactions
intersect(bioplanet[["EGFR1 pathway"]], bioplanet[["Interleukin-2 signaling pathway"]]) 
# 36 genes in common, 0.03738318 Weight.bp
intersect(bioplanet[["Axon guidance"]], bioplanet[["Ephrin receptor B forward pathway"]]) 
# 32 genes in common, 0.09065156 Weight.bp
intersect(bioplanet[["Biosynthesis of unsaturated fatty acids"]], bioplanet[["Lysosome"]]) 
# 0, 0
# Attempt to plot this - 
thresh.vec2 <- c(seq(0.2, 0.082, -0.002), seq(0.08, 0.01, -0.01))
identical (thresh.vec2, net.atts.df$threshold) # TRUE
dev.new()
plot(pathway.net.list4[[68]]$Weight.clust ~ pathway.net.list4[[68]]$Weight.bp, pch=16, col=alpha("grey88", 0.2), xlab="bioplanet weight", ylab="PTM cluster weight")
plotcolors <- rainbow(68)
for (i in 1:length(pathway.net.list4)){
  points(head(pathway.net.list4[[i]]$Weight.clust, 20) ~ head(pathway.net.list4[[i]]$Weight.bp, 20), pch=16, col=alpha(plotcolors[i], 0.5))
  print(i)
  Sys.sleep(0.04)
  }
# This is interesting, the top clusters are the same for the first set, and have lower PTM cluster weight. 
# At high threshold, a few clusters contribute to pathway-pathway edge weight. 
# At lower threshold, more clusters contribute, the total weight is increased, and different Edges are selected. 
selected.pathway.edges <- big.path.net[big.path.net$Weight.clust > 0.8 & big.path.net$Weight.bp < 0.001,]
# 251 pathway interactions at this setting
# Graph
spe.nodes <- data.frame(id=unique(c(selected.pathway.edges$source, selected.pathway.edges$target)))
look.suid <- createNetworkFromDataFrames(spe.nodes, selected.pathway.edges[ c("source", "target", "Weight", "interaction")], title="Pathway Interactions Weight.clust > 0.8 & Weight.bp < 0.001", collection = "Interactions")
layoutNetwork("force-directed")
# Tweak settings
selected.pathway.edges2 <- big.path.net[big.path.net$Weight.clust > 1.6 & big.path.net$Weight.bp < 0.01,]
# 126 pathway interactions at this setting
# Graph
spe.nodes2 <- data.frame(id=unique(c(selected.pathway.edges2$source, selected.pathway.edges2$target)))
look2.suid <- createNetworkFromDataFrames(spe.nodes2, selected.pathway.edges2[ c("source", "target", "Weight", "interaction")], title="Pathway Interactions, Weight > 1.6; Weight.bp < 0.01", collection = "Interactions")
layoutNetwork("genemania-force-directed")
intersect(bioplanet[["Protein processing in the endoplasmic reticulum"]], bioplanet[["Developmental biology"]]) 
# [1] "HSP90AA1" "HSP90AB1"  **!
#
# Tweak settings
selected.pathway.edges3 <- big.path.net[big.path.net$Weight.clust > 1.1 & big.path.net$Weight.bp < 0.01,]
# 318 pathway interactions at this setting
# Graph
spe.nodes3 <- data.frame(id=unique(c(selected.pathway.edges3$source, selected.pathway.edges3$target)))
look3.suid <- createNetworkFromDataFrames(spe.nodes3, selected.pathway.edges3[ c("source", "target", "Weight", "interaction")], title="Pathway Interactions, Weight > 1.1; Weight.bp < 0.01", collection = "Interactions")
layoutNetwork("genemania-force-directed")
print(unlist(spe.nodes3$id))
# Think about intersection with clusters from:
selected.pathway.edges4 <- big.path.net[big.path.net$Weight.clust > 1 & big.path.net$Weight.bp < 1,]
# 1045 pathway interactions at this setting
# Graph
spe.nodes4 <- data.frame(id=unique(c(selected.pathway.edges4$source, selected.pathway.edges4$target)))
spe.nodes4 <- data.frame(id=unique(c(selected.pathway.edges4$source, selected.pathway.edges4$target)))
look4.suid <- createNetworkFromDataFrames(spe.nodes4, selected.pathway.edges4[ c("source", "target", "Weight", "interaction")], title="Pathway Interactions, Weight > 1; Weight.bp < 1", collection = "Interactions")
layoutNetwork("genemania-force-directed")
# EAxon Guidance Pathway
# zymes.bioplanet.common, zymes.gzclusters.bioplanet.common, zymes.gzclusters.common, zymes.list 
#

# Explore example from BioPlanetTest for links between EGFR and FASN
a <- "Signaling by EGFR in cancer"
b <- "Fatty acid, triacylglycerol, and ketone body metabolism"
filter.edges.0(c(a,b), big.path.net)[,c(3,6,7)]
# Weight.bp   Weight.clust     Weight
# 0.003533569 0.02083333      0.01729976
# Not high.
# Where else is FASN?
fasn.bp <- bioplanet[grep("FASN", bioplanet)]
fasn.paths <- names (fasn.bp)
egfr.bp <- bioplanet[grep("EGFR", bioplanet)]
intersect(unlist(fasn.paths), unlist(egfr.paths))
egfr.paths <- names (egfr.bp)
# Examine interactions
egfr.fasn.between <- filter.edges.between(egfr.paths, fasn.paths, big.path.net)
# Gets all edges connecting these: 410
egfr.fasn <- filter.edges.0(c(egfr.paths, fasn.paths), big.path.net) 
# 1463 edges, gets everything in between all nodes.
# Examine weights
egfr.fasn.between[egfr.fasn.between$Weight>1,]
ef.sel <- egfr.fasn.between[egfr.fasn.between$Weight.clust > 0.5 & egfr.fasn.between$Weight.bp < 0.01,]
# "Metabolism" is huge! > 1600 genes
efnomet <- egfr.fasn.between[-grep("Metabolism",  egfr.fasn.between$source),]
efnomet <- efnomet[-grep("Metabolism",  efnomet$target),]
# 335 left
efnomet[efnomet$Weight.clust > 0.3 & efnomet$Weight.bp < 0.01,]
# Disease is too big
efnomet <- efnomet[-grep("Disease",  efnomet$source),]
efnomet <- efnomet[-grep("Disease",  efnomet$target),]
# 321 left *** looking more interesting
hist(efnomet$Weight, col="blue", , breaks=50)
efselect <- efnomet[efnomet$Weight>=0.2,] # 37 edges
efselect.nodes <- data.frame(id=unique(c(efselect$source, efselect$target)))
look3.suid <- createNetworkFromDataFrames(efselect.nodes, efselect[ c("source", "target", "Weight", "interaction")], title="EGFR FASN Pathway Interactions", collection = "Interactions")
layoutNetwork("genemania-force-directed")
setEdgeWidths.RCy32(efselect, factor=5, log=F)
#
# Other Pathway interactions 
egfnames <- names(bioplanet)[grep("EGF", names(bioplanet))] %w/o% names(bioplanet)[grep("VEGF", names(bioplanet))]
names(bioplanet)[grep("Endocytosis", names(bioplanet))]
filter.edges.between("Endocytosis", egfnames, selected.pathway.edges)
#####>>>>-----
# Plots for figures
# Examine above pathway genes as in Karenguolin12
a <- bioplanet[["Protein processing in the endoplasmic reticulum"]]
b <- bioplanet[["Developmental biology"]]
ab <- intersect (a, b)
ab.all <- unique(c(a, b)) # 584
# only "HSP90AA1" "HSP90AB1" 
gzpathgenes <- ab.all[ab.all %in% gzallt.gene.key$Gene.Name] # 201
ER_Dev <- composite.shortest.paths(genes1=a[a %in% gzallt.gene.key$Gene.Name], genes2=b[b %in% gzallt.gene.key$Gene.Name], network=gzalltgene.physical.cfn.merged, exclude="MYH9")
# 
# 2520  merged edges
test.ld <- composite.shortest.paths(genes1=a[a %in% ld.gene.key$Gene.Name], genes2=b[b %in% ld.gene.key$Gene.Name], network=ldgene.physical.cfn.merged, exclude="")
# 1320 edges
# new graphs
# graph.cfn.cccn (test.gz, ld=FALSE, gz=TRUE, only.cfn=TRUE)
cccn1 <- graph.cfn.cccn (ER_Dev, ld=FALSE, gz=TRUE, only.cfn=FALSE)
all.ratio.styles()
toggleGraphicsDetails()
# Follow up top EGFR FASN interaction from efselect above
a <- bioplanet[["Lipid and lipoprotein metabolism"]]
b <- bioplanet[["Developmental biology"]]
ab <- intersect (a, b) # 58 - loads of MED mediator proteins!
b <- bioplanet[["Lipid and lipoprotein metabolism"]]
ab <- intersect (a, b) # 489!
# Look for somethihng from Pathway Interactions Weight.clust > 0.8 & Weight.bp < 0.001 selected.pathway.edges
filter.edges.1("EGF/EGFR signaling pathway", selected.pathway.edges)
# None have any bp interactions. Choose Transmembrane transport of small molecules (druggable?)
a <- bioplanet[["Transmembrane transport of small molecules"]]
b <- bioplanet[["EGF/EGFR signaling pathway"]]
ab <- intersect (a, b) #0
ab.all <- unique(c(a, b)) # 573
# 
intersect(bioplanet[["Protein processing in the endoplasmic reticulum"]], bioplanet[["EGF/EGFR signaling pathway"]]) #[1] "MAPK8" "MAPK9"
intersect(bioplanet[["Apoptosis regulation"]], bioplanet[["EGF/EGFR signaling pathway"]])
# none
gzpathgenes <- ab.all[ab.all %in% gzallt.gene.key$Gene.Name] # 201
egfr.transporters <- composite.shortest.paths(genes1=a[a %in% gzallt.gene.key$Gene.Name], genes2=b[b %in% gzallt.gene.key$Gene.Name], network=gzalltgene.physical.cfn.merged, exclude="MYH9")
dim(egfr.transporters) # 2115 4
cccn2 <- graph.cfn.cccn (egfr.transporters, ld=FALSE, gz=TRUE, only.cfn=FALSE)
cfn2 <- graph.cfn.cccn (egfr.transporters, ld=FALSE, gz=TRUE, only.cfn=TRUE)
# Still very complex
FixEdgeDprops.RCy32()
all.ratio.styles()
toggleGraphicsDetails()
# Save as Transmembrane_EGF_CFN_CCCN.cys
# What are the edges *between* uniuqe pathway genes?
look4 <- filter.edges.between( bioplanet[["Transmembrane transport of small molecules"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
# 13 edges
selectNodes(nodes=extract.gene.names.RCy3(look4), by.col="id")
selectFirstNeighbors()
createSubnetwork(subnetwork.name="Edges Between Pathays")
# Still to complex, try this (again)
cccn4 <- graph.cfn.cccn (look4, ld=FALSE, gz=TRUE, only.cfn=FALSE)
# Try mulitple times, then
delete.bad.networks()
# Experimental function to check graphs as you go, up to 10 tries
cccn4 <- graph.cfn.cccn.check (look4, ld=FALSE, gz=TRUE, only.cfn=FALSE)

#setCorrEdgeAppearance(cccn4) 
#setCorrEdgeAppearance(getTableColumns('edge'))
#FixEdgeDprops.RCy32()
look5 <- filter.edges.between( bioplanet[["Protein metabolism"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
# 32 edges
get.all.gene.names.from.peps(unique(c(look5$source, look5$target)))
cccn5 <- graph.cfn.cccn (look5, ld=FALSE, gz=TRUE, only.cfn=FALSE)
setCorrEdgeAppearance(cccn5) 
setCorrEdgeAppearance(getTableColumns('edge'))
FixEdgeDprops.RCy32()
# Merge these two for styles
look45 <- rbind(look4, look5)
# merged these two networks manualy in Cytoscape
nodeDprops.RCy32()
all.ratio.styles()
# ----This one has Weight.Clust 4.5 and Weight.bp 0.007
filter.edges.between("Metabolism", "EGF/EGFR signaling pathway", total.pathway.net)
intersect(bioplanet[["Metabolism"]], bioplanet[["EGF/EGFR signaling pathway"]])
# 13 genes including PI3 kinases and PI phosphatases, CAV1, AKT1, PLCG1, PLD1, IQGAP1,...
look6 <- filter.edges.between( bioplanet[["Metabolism"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
# 68 edges
get.all.gene.names.from.peps(unique(c(look6$source, look6$target)))
cccn6 <- graph.cfn.cccn (look6, ld=FALSE, gz=TRUE, only.cfn=FALSE)

# Axon guidance is related to metastasis
intersect(bioplanet[["Axon guidance"]], bioplanet[["EGF/EGFR signaling pathway"]])
look7 <- filter.edges.between( bioplanet[["Axon guidance"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
# 168 edges
get.all.gene.names.from.peps(unique(c(look6$source, look6$target)))
cccn7 <- graph.cfn.cccn (look7, ld=FALSE, gz=TRUE, only.cfn=FALSE)
#

# Try mulitple times, then
delete.bad.networks()
# Merge all cccns and delete styles, then redo them
mystyles <- getVisualStyleNames()
deleteVisualStyle(mystyles[grep("tyle", mystyles)])
nodeDprops.RCy32()
all.ratio.styles()
setCorrEdgeAppearance(getTableColumns('edge'))
FixEdgeDprops.RCy32()

# Who is in which pathway?
look4.df <- data.frame(Gene.Name=sort(unique(c(look4$source, look4$target))))
look4.df$EGFR.path <- sapply(look4.df$Gene.Name, function (x) x %in% b)
look4.df$Transporter.path <- sapply(look4.df$Gene.Name, function (x) x %in% a)
#__________________
gzpathgenes <- ab.all[ab.all %in% gzallt.gene.key$Gene.Name] # 201
ER_Dev <- composite.shortest.paths(genes1=a[a %in% gzallt.gene.key$Gene.Name], genes2=b[b %in% gzallt.gene.key$Gene.Name], network=gzalltgene.physical.cfn.merged, exclude="MYH9")
# 
# 2520  merged edges
test.ld <- composite.shortest.paths(genes1=a[a %in% ld.gene.key$Gene.Name], genes2=b[b %in% ld.gene.key$Gene.Name], network=ldgene.physical.cfn.merged, exclude="")
# 1320 edges
# new graphs
# graph.cfn.cccn (test.gz, ld=FALSE, gz=TRUE, only.cfn=TRUE)
cccn1 <- graph.cfn.cccn (ER_Dev, ld=FALSE, gz=TRUE, only.cfn=FALSE)
all.ratio.styles()
toggleGraphicsDetails()
# This is useful to see which clusters are predominantly affected by particular drugs. 
# What are the edges *between* uniuqe pathway genes?
ernodev <- a %w/o% b
devnoer <- b %w/o% a
look5 <- filter.edges.between(ernodev, devnoer, edge.file=gzalltgene.physical.cfn.merged)
# 34 edges, including CLTC, EGFR, MAPK9 (JNK2)...
cccn5 <- graph.cfn.cccn (look5, ld=FALSE, gz=TRUE, only.cfn=FALSE)
# Note: run all.ratio.styles() if not run above.
# Edges are messed up: 'shared interaction' is wrong!
# Doesnt help: 
edgeDprops.RCy32(); all.ratio.styles() 
# Replotting also doesn't help. Try restarting cytoscape: This worked for some but not all edges. Glitch!
# Cytoscape 3.8 doesn't start properly from .sh script. Open manually, try again
all.ratio.styles()
toggleGraphicsDetails()

#-------------------------------------------
# A cloud
plot(total.pathway.net$Combined.Weight~total.pathway.net$Weight.clust, log="xy", pch=20, col=alpha(col2hex("magenta"), 0.1), cex=1.8)
points(big.path.net$Weight~big.path.net$Weight.clust,  pch=19, col=alpha(col2hex("blue"), 0.1), cex=1.8)
points(total.pathway.net$Weight.normalized~total.pathway.net$Weight.clust,  pch=20, col=alpha(col2hex("yellow"), 0.1), cex=1.8)
# Try different order
dev.new()
plot(total.pathway.net$Weight.normalized~total.pathway.net$Weight.clust, ylab= "Transformed weight", xlab="Cluster evidence weight", pch=20, col=alpha(col2hex("green"), 0.1), log="xy", cex=1.8)
points(total.pathway.net$Combined.Weight~total.pathway.net$Weight.clust,  pch=20, col=alpha(col2hex("magenta"), 0.1), cex=1.8)
points(big.path.net$Weight~big.path.net$Weight.clust,  pch=20, col=alpha(col2hex("blue"), 0.1), cex=1.8)
# ***
legend("bottomright", pt.cex=1.8, pch=20, col=c("green", "magenta", "blue"), legend = c("Weight normalized", "Weight combined", "Weight normalized, filtered"))
# Without log
dev.new()
plot(total.pathway.net$Weight.normalized~total.pathway.net$Weight.clust, ylab= "Transformed weight", xlab="Cluster evidence weight", pch=20, col=alpha(col2hex("green"), 0.25), cex=1.8)
points(total.pathway.net$Combined.Weight~total.pathway.net$Weight.clust,  pch=20, col=alpha(col2hex("magenta"), 0.2), cex=1.8)
points(big.path.net$Weight~big.path.net$Weight.clust,  pch=20, col=alpha(col2hex("blue"), 0.2), cex=1.8)
# ***
legend("bottomright", pt.cex=1.8, pch=20, col=c("green", "magenta", "blue"), legend = c("Weight normalized", "Weight combined", "Weight normalized, filtered"))
# ***
hist (total.pathway.net$Weight.clust, col="green", breaks=40000, xlim=c(0, 0.05))
hist (big.path.net$Weight.clust, col="red", breaks=40000, xlim=c(0, 0.05), add=TRUE)
# ***
dev.new()
plot(density(total.pathway.net$Weight.clust), col="green3", xlim=c(0,2))
lines(density(big.path.net$Weight.clust), col="red",  xlim=c(0,2))
lines(density(total.pathway.net$Combined.Weight), col="blue",  xlim=c(0, 2))

# zoom in:
plot(density(total.pathway.net$Weight.clust), col="green3", xlim=c(0,2), ylim=c(0, 1))
lines(density(big.path.net$Weight.clust), col="red",  xlim=c(0,2), ylim=c(0, 1))
plot(density(total.pathway.net$Weight.clust), col="green3", xlim=c(2,5), ylim=c(0, 0.02))
lines(density(big.path.net$Weight.clust), col="red",  xlim=c(2,5), ylim=c(0, 0.02))
lines(density(total.pathway.net$Combined.Weight), col="blue",  xlim=c(2,5), ylim=c(0, 0.02))
# What does the graph look like? Repeat the above
# Focus on cluster evidence but also graph bp edges
tpn <- total.pathway.net.no.bp
tpn.small <- tpn[tpn$Weight>2,]
# 249 edges
# Plots
hist(tpn$Weight, col='orange', breaks=50) # Max 6.3
# Now get the other edges
tpn.small.nodes <- data.frame(id=unique(c(tpn.small$source, tpn.small$target)))
# 124 nodes
pcn.small <- filter.edges.0(tpn.small.nodes$id, pathway.crosstalk.network)
# This gives all edges: 12116. Need to filter.
bp.small <- filter.edges.0(tpn.small.nodes$id, bioplanetjaccardedges)
# 4490 edges - need to filter
hist(bp.small$Weight, breaks=100, col="yellow")
bp.small.1 <- bp.small[bp.small$Weight>0.2,] # 280 edges
bp.small.2 <- bp.small[bp.small$Weight>0.3,] # 186 edges
bp.small.3 <- bp.small[bp.small$Weight>0.4,] # 113 edges
# Bind together
pcn.small <- rbind(tpn.small[,c("source", "target", "Weight", "interaction")], bp.small.2)
# Graph
look8.suid <- createNetworkFromDataFrames.check(tpn.small.nodes, pcn.small, title="Pathway Interactions, Weight > 2.0; Weight.bp graphed", collection = "Interactions")
layoutNetwork("genemania-force-directed")
setEdgeSelectionColorDefault (col2hex("chartreuse"))
edgeColors <- c(col2hex("purple"), col2hex("green"))
edgeTypes <- c("cluster evidence", "pathway Jaccard similarity")
setEdgeColorMapping( 'interaction', edgeTypes, edgeColors, 'd', default.color="#FFFFFF")
style.name <- "PCN style 3"
copyVisualStyle('default', style.name)
setVisualStyle(style.name)
# Get a subset of nodes using Cytoscape
funnodes <- getAllNodes()
# Metabolism, Axon guidance, Protein metabolism, Messenger RNA processing, EGF/EGFR signaling pathway, Apoptosis regulation, Translation, Spliceosome, Protein processing in the endoplasmic reticulum, Transmembrane transport of small molecules, SLC-mediated transmembrane transport
# With new weights, the transporter-EGF signaling weight is different
filter.edges.between("Transmembrane transport of small molecules", "EGF/EGFR signaling pathway", pathway.crosstalk.network)
# Weight 1.66
filter.edges.between("Endocytosis", "EGF/EGFR signaling pathway", pathway.crosstalk.network)
# Weight cluster 2.09728951; bp 0.086
# Try picking some from this list and filtering edges
selections <- funnodes[c(2,3,6,7,12,13,15,17,19,20,21,23,27,34,46,49,54,56,58)]
selections.nodes <- data.frame(id=selections)
selection.edges <- filter.edges.0(selections, pathway.crosstalk.network)
# 269 edges
hist(selection.edges$Weight, breaks=50, col="violet")
hist(selection.edges$Weight, breaks=100, col="darkviolet", xlim=c(0,1))
look9.suid <- createNetworkFromDataFrames.check(selections.nodes, selection.edges, title="Selected Pathway Interactions, unfiltered", collection = "Interactions")
layoutNetwork("genemania-force-directed")
# Too many edges!
# Filter 
selected.edges.bp <- selection.edges[selection.edges$interaction=="pathway Jaccard similarity",]
selected.edges.cpe <- selection.edges[selection.edges$interaction=="cluster evidence",]
hist(selected.edges.bp$Weight, breaks=100, col="darkviolet", xlim=c(0,0.1))
hist(selected.edges.cpe$Weight, breaks=100, col="darkviolet", xlim=c(0,2))
sel.cpe.filtered <- selected.edges.cpe[selected.edges.cpe$Weight>1,] # 141
sel.bp.filtered <- selected.edges.bp[selected.edges.bp$Weight>0.01,] # 59
sel.filtered <- rbind(sel.cpe.filtered, sel.bp.filtered)
sel.nodes <- data.frame(id=unique(c(sel.filtered$source, sel.filtered$target)))
look10.suid <- createNetworkFromDataFrames.check(sel.nodes, sel.filtered, title="Selected Pathway Interactions, clust weight>1, bp weight>0.01", collection = "Interactions")
layoutNetwork("genemania-force-directed")
# Still too many edges
sel.cpe.filtered2 <- selected.edges.cpe[selected.edges.cpe$Weight>2,] # 83
sel.bp.filtered2 <- selected.edges.bp[selected.edges.bp$Weight>0.035,] # 36
sel.filtered2 <- rbind(sel.cpe.filtered2, sel.bp.filtered2)
sel.nodes2 <- data.frame(id=unique(c(sel.filtered2$source, sel.filtered2$target)))
look10.suid <- createNetworkFromDataFrames.check(sel.nodes2, sel.filtered2, title="Selected Pathway Interactions, clust weight>2, bp weight>0.035", collection = "Interactions")
layoutNetwork("genemania-force-directed")
# Filter further using no.bp
sel.cpe.filtered3 <- filter.edges.0(selections, total.pathway.net.no.bp) # 73
sel.cpe.filtered3 <-  sel.cpe.filtered3[sel.cpe.filtered3$Weight>1,] # 56
sel.bp.filtered3 <- selected.edges.bp[selected.edges.bp$Weight>0.039,] # 34
sel.bp.filtered3a <- sel.cpe.filtered3[,c("source", "target", "Weight", "interaction")]
sel.bp.filtered3a$interaction <- "cluster evidence"
sel.filtered3 <- rbind(sel.bp.filtered3a, sel.bp.filtered3)
sel.nodes3 <- data.frame(id=unique(c(sel.filtered3$source, sel.filtered3$target)))
look12.suid <- createNetworkFromDataFrames.check(sel.nodes3, sel.filtered3, title="Selected Pathway Interactions, clust weight>1, bp weight>0.039", collection = "Interactions")
layoutNetwork("genemania-force-directed")
# *** starting to look managable.
#****
# Explore example from BioPlanetTest for links between EGFR and FASN
a <- "Signaling by EGFR in cancer"
b <- "Fatty acid, triacylglycerol, and ketone body metabolism"
filter.edges.0(c(a,b), big.path.net)[,c(3,6,7)]
# Weight.bp   Weight.clust     Weight
# 0.003533569 0.02083333      0.01729976
# Not high.
# Where else is FASN?
fasn.bp <- bioplanet[grep("FASN", bioplanet)]
fasn.paths <- names (fasn.bp)
egfr.bp <- bioplanet[grep("EGFR", bioplanet)]
intersect(unlist(fasn.paths), unlist(egfr.paths))
egfr.paths <- names (egfr.bp)
# Examine interactions
egfr.fasn.between <- filter.edges.between(egfr.paths, fasn.paths, big.path.net)
# Gets all edges connecting these: 410
egfr.fasn <- filter.edges.0(c(egfr.paths, fasn.paths), big.path.net) 
# 1463 edges, gets everything in between all nodes.
# Examine weights
egfr.fasn.between[egfr.fasn.between$Weight>1,]
ef.sel <- egfr.fasn.between[egfr.fasn.between$Weight.clust > 0.5 & egfr.fasn.between$Weight.bp < 0.01,]
# "Metabolism" is huge! > 1600 genes
efnomet <- egfr.fasn.between[-grep("Metabolism",  egfr.fasn.between$source),]
efnomet <- efnomet[-grep("Metabolism",  efnomet$target),]
# 335 left
efnomet[efnomet$Weight.clust > 0.3 & efnomet$Weight.bp < 0.01,]
# Disease is too big
efnomet <- efnomet[-grep("Disease",  efnomet$source),]
efnomet <- efnomet[-grep("Disease",  efnomet$target),]
# 321 left *** looking more interesting
hist(efnomet$Weight, col="blue", , breaks=50)
efselect <- efnomet[efnomet$Weight>=0.2,] # 37 edges
efselect.nodes <- data.frame(id=unique(c(efselect$source, efselect$target)))
look3.suid <- createNetworkFromDataFrames(efselect.nodes, efselect[ c("source", "target", "Weight", "interaction")], title="EGFR FASN Pathway Interactions", collection = "Interactions")
layoutNetwork("genemania-force-directed")
setEdgeWidths.RCy32(efselect, factor=5, log=F)
#
# Other Pathway interactions 
egfnames <- names(bioplanet)[grep("EGF", names(bioplanet))] %w/o% names(bioplanet)[grep("VEGF", names(bioplanet))]
names(bioplanet)[grep("Endocytosis", names(bioplanet))]
filter.edges.between("Endocytosis", egfnames, selected.pathway.edges)
#####>>>>-----
# Plots for figures
# Examine above pathway genes as in Karenguolin12
a <- bioplanet[["Protein processing in the endoplasmic reticulum"]]
b <- bioplanet[["Developmental biology"]]
ab <- intersect (a, b)
ab.all <- unique(c(a, b)) # 584
# only "HSP90AA1" "HSP90AB1" 
gzpathgenes <- ab.all[ab.all %in% gzallt.gene.key$Gene.Name] # 201
ER_Dev <- composite.shortest.paths(genes1=a[a %in% gzallt.gene.key$Gene.Name], genes2=b[b %in% gzallt.gene.key$Gene.Name], network=gzalltgene.physical.cfn.merged, exclude="MYH9")
# 
# 2520  merged edges
test.ld <- composite.shortest.paths(genes1=a[a %in% ld.gene.key$Gene.Name], genes2=b[b %in% ld.gene.key$Gene.Name], network=ldgene.physical.cfn.merged, exclude="")
# 1320 edges
# new graphs
# graph.cfn.cccn (test.gz, ld=FALSE, gz=TRUE, only.cfn=TRUE)
cccn1 <- graph.cfn.cccn (ER_Dev, ld=FALSE, gz=TRUE, only.cfn=FALSE)
all.ratio.styles()
toggleGraphicsDetails()
# Follow up top EGFR FASN interaction from efselect above
a <- bioplanet[["Lipid and lipoprotein metabolism"]]
b <- bioplanet[["Developmental biology"]]
ab <- intersect (a, b) # 58 - loads of MED mediator proteins!
b <- bioplanet[["Lipid and lipoprotein metabolism"]]
ab <- intersect (a, b) # 489!
# Look for something from Pathway Interactions Weight.clust > 0.8 & Weight.bp < 0.001 selected.pathway.edges
filter.edges.1("EGF/EGFR signaling pathway", selected.pathway.edges)
# None have any bp interactions. Choose Transmembrane transport of small molecules (druggable?)
a <- bioplanet[["Transmembrane transport of small molecules"]]
b <- bioplanet[["EGF/EGFR signaling pathway"]]
ab <- intersect (a, b) #0
ab.all <- unique(c(a, b)) # 573
# 
gzpathgenes <- ab.all[ab.all %in% gzallt.gene.key$Gene.Name] # 201
egfr.transporters <- composite.shortest.paths(genes1=a[a %in% gzallt.gene.key$Gene.Name], genes2=b[b %in% gzallt.gene.key$Gene.Name], network=gzalltgene.physical.cfn.merged, exclude="MYH9")
cccn2 <- graph.cfn.cccn.check (egfr.transporters, ld=FALSE, gz=TRUE, only.cfn=FALSE)
#FixEdgeDprops.RCy32()
all.ratio.styles()
toggleGraphicsDetails()
# What are the edges *between* uniuqe pathway genes?
look4 <- filter.edges.between( bioplanet[["Transmembrane transport of small molecules"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
# 13 edges
selectNodes(nodes=extract.gene.names.RCy3(look4), by.col="id")
selectFirstNeighbors()
createSubnetwork(subnetwork.name="Edges Between Pathways")
# Still to complex, try this (again)
cccn4 <- graph.cfn.cccn (look4, ld=FALSE, gz=TRUE, only.cfn=FALSE)
setCorrEdgeAppearance(cccn4) 
setCorrEdgeAppearance(getTableColumns('edge'))
FixEdgeDprops.RCy32()
# Already done all.ratio.styles()
# **** Got two network figures before I couldn't fix edges! Try manually:
# Who is in which pathway?
look4.df <- data.frame(Gene.Name=sort(unique(c(look4$source, look4$target))))
look4.df$EGFR.path <- sapply(look4.df$Gene.Name, function (x) x %in% b)
look4.df$Transporter.path <- sapply(look4.df$Gene.Name, function (x) x %in% a)
# >>>>>>>>>>>>>>>>>
# "Friends of friends" approach
# Use CFN here; keep genetic interactions (gzalltgene.cfn.rcy3)
a <- bioplanet[["Transmembrane transport of small molecules"]]
b <- bioplanet[["EGF/EGFR signaling pathway"]]
transporter.friends <- filter.edges.1(a, gzalltgene.cfn.rcy3) 
EGFR.friends <- filter.edges.1(b, gzalltgene.cfn.rcy3) 
friends.of.friends <- intersect(extract.gene.names(transporter.friends), extract.gene.names(EGFR.friends)) %w/o% c(a,b)
# 120 mutual friends
# Network of these
fof.cfn <- filter.edges.0(c(a,b,friends.of.friends),gzalltgene.cfn.rcy3)
# 2634 edges
cfn1 <- graph.cfn.cccn (fof.cfn, ld=FALSE, gz=TRUE, only.cfn=TRUE)
toggleGraphicsDetails()
# Use simpler merged version and filter.edges.between
fof.cfn2 <- filter.edges.between(c(a,b),friends.of.friends, gzalltgene.physical.cfn.merged) # now 296 edges
cfn2 <- graph.cfn.cccn (fof.cfn2, ld=FALSE, gz=TRUE, only.cfn=TRUE)
# separate out
selectNodes(nodes=a, by.col="id", preserve.current.selection = F)
selectNodes(nodes=b, by.col="id", preserve.current.selection = F)
# ****
# Try this with different views, including composite shortest paths and physical between in "EGF transporters.cys"
selectNodes(nodes=friends.of.friends, by.col="id", preserve.current.selection = F)
# Who is not in friends.of.friends in composite shortest paths: egfr.transporters
csp.genes <- outersect(extract.gene.names(egfr.transporters), c(a, b, friends.of.friends)) # 253 genes§
selectNodes(nodes=csp.genes, by.col="id", preserve.current.selection = F)
# look at CFN version
selectNodes(nodes=extract.gene.names(egfr.transporters), by.col="id", preserve.current.selection = F)
createSubnetwork(subnetwork.name="EGFR transporters composite shortest paths CFN")
# Back to csp for CCCN version
selectNodes(nodes=extract.gene.names(egfr.transporters), by.col="id", preserve.current.selection = F)
invertNodeSelection()
createSubnetwork(subnetwork.name="EGFR transporters composite shortest paths CCCN")
# Idea: betweenness of nodes between pathways
egfr.transporters.g <- graph_from_data_frame(egfr.transporters)
e.t.between <- betweenness(egfr.transporters.g)
e.t.b.df <- data.frame(Gene.Name=names(e.t.between), betweenness=e.t.between)
e.t.b.df <- e.t.b.df[order(e.t.b.df$betweenness, decreasing=TRUE),]
# ***
# Who is where? And in what family
between.genes <- e.t.b.df[e.t.b.df$betweenness>0, "Gene.Name"]
intersect(between.genes,   bioplanet[["Transmembrane transport of small molecules"]]) # 12
intersect(between.genes,   bioplanet[["EGF/EGFR signaling pathway"]]) # 12
# 37; EGFR is at the top. Are these related to overall betweenness?
etb.df <- left_join(e.t.b.df, gzgene.cfn.netatts, by="Gene.Name")
gz.cf.ppibetween <- gz.cf[which(gz.cf$Node.ID=="gene"), c("Gene.Name", "ppibetween", "norm.ppibetween", "HPRD.Function")]
etb.df <- left_join(etb.df, gz.cf.ppibetween[, c("Gene.Name", "HPRD.Function")], by="Gene.Name")
etb.df$enzyme.class <- sapply(etb.df$Gene.Name, function (x) names(zymes.list[grep(x, zymes.list)]))

# Relationships to pathway crosstalk network betweenness from CFN and PPI networks
plot(etb.df$betweenness~etb.df$ppibetween, pch=19, col="seagreen")
summary(lm(etb.df$betweenness~etb.df$ppibetween))
# R-squared:  0.2305 
plot(etb.df$betweenness~etb.df$norm.ppibetween, pch=19, col="forestgreen")
# norm.ppibetween is CFN ppibetween/allppibetween
summary(lm(etb.df$betweenness~etb.df$norm.ppibetween))
# R-squared: 0.02623
plot(etb.df$betweenness~etb.df$allppibetween, pch=19, col="forestgreen")
summary(lm(etb.df$betweenness~etb.df$allppibetween))
# R-squared: 0.07327

# Which genes are high betweenness relative to total CFN?
#
etb.df[which(etb.df$betweenness>etb.df$allppibetween), "Gene.Name"]
# [1] "CTTN"  "KRT18"
etb.df[which(etb.df$betweenness>etb.df$ppibetween), "Gene.Name"] -> highbetweeners
#[1] "CTTN"   "EFNB1"  "OCLN"   "PAK2"   "GJA1"   "BCAR1"  "PRKACA" "PIK3R1" "MAPK14" "EPHB4" "RALA"   "PLAA"   "IL6ST"  "PSAT1"  "KIT"   

selectNodes(highbetweeners, by.col="id", preserve.current.selection = F)
selectFirstNeighbors()
# Ask if these are enzymes
etb.df[etb.df$Gene.Name %in% highbetweeners, c("Gene.Name", "HPRD.Function", "enzyme.class")]
# Group by enzyme or function
etb.RNA <- etb.df[grep("RNA", etb.df[,c("HPRD.Function","enzyme.class")], fixed = T), "Gene.Name"]
etb.kinase <- etb.df[grep("kinase", etb.df[,c("HPRD.Function","enzyme.class")]), "Gene.Name"]

# For figures with no node labels
setNodeLabelColorDefault(col2hex("transparent"), style.name = "H1781_AfatinibRatio Style",)

#__________________
gzpathgenes <- ab.all[ab.all %in% gzallt.gene.key$Gene.Name] # 201
ER_Dev <- composite.shortest.paths(genes1=a[a %in% gzallt.gene.key$Gene.Name], genes2=b[b %in% gzallt.gene.key$Gene.Name], network=gzalltgene.physical.cfn.merged, exclude="MYH9")
# 
# 2520  merged edges
test.ld <- composite.shortest.paths(genes1=a[a %in% ld.gene.key$Gene.Name], genes2=b[b %in% ld.gene.key$Gene.Name], network=ldgene.physical.cfn.merged, exclude="")
# 1320 edges
# new graphs
# graph.cfn.cccn (test.gz, ld=FALSE, gz=TRUE, only.cfn=TRUE)
cccn1 <- graph.cfn.cccn (ER_Dev, ld=FALSE, gz=TRUE, only.cfn=FALSE)
all.ratio.styles()
toggleGraphicsDetails()
# This is useful to see which clusters are predominantly affected by particular drugs. 
# What are the edges *between* uniuqe pathway genes?
ernodev <- a %w/o% b
devnoer <- b %w/o% a
look5 <- filter.edges.between(ernodev, devnoer, edge.file=gzalltgene.physical.cfn.merged)
# 34 edges, including CLTC, EGFR, MAPK9 (JNK2)...
cccn5 <- graph.cfn.cccn (look5, ld=FALSE, gz=TRUE, only.cfn=FALSE)
# Note: run all.ratio.styles() if not run above.
# Edges are messed up: 'shared interaction' is wrong!
# Doesnt help: 
edgeDprops.RCy32(); all.ratio.styles() 
# Replotting also doesn't help. Try restarting cytoscape: This worked for some but not all edges. Glitch!
# Cytoscape 3.8 doesn't start properly from .sh script. Open manually, try again
all.ratio.styles()
toggleGraphicsDetails()

#-------------------------------------------
dev.new()
plot(unlist(density.list2) ~ thresh.vec2, col="gray100", pch=16, xlab="Cluster Evidence Threshold", ylab="Pathway Network Density", cex.lab=1.25)
net.nodes.index <- (net.nodes/10000)
net.edges.index <- (net.edges/1000000)
#symbols(unlist(density.list) ~ thresh.vec,  circles=log(net.nodes)/500, inches = FALSE, bg="deepskyblue", fg = "darkblue",  add=TRUE)
# symbols(unlist(density.list) ~ thresh.vec,  circles=log(net.edges)/1000, inches = FALSE, bg="orange", fg = "red",  add=TRUE)
symbols(unlist(density.list2) ~ thresh.vec2,  circles=net.nodes.index/5, inches = FALSE, bg="deepskyblue", fg = "darkblue",  add=TRUE)

text(topallppi[, "allppibetween"], topallppi[, "ppibetween"], col="red", cex=0.8, labels=noquote(as.character(topallppi $Gene.Name)))
text(topesslhubs[, "allppibetween"], topesslhubs[, "ppibetween"], col="black", cex=0.8, labels=noquote(as.character(topesslhubs$Gene.Name)))  

bioplanetjaccard.g <- bioplanet.g

# Normalize: put two weights on same scale to make Weight.clust from 0 to 1
tpnn <- total.pathway.net[,c(1:4,6,5)]
tpnn$Weight.clust <- tpnn$Weight.clust/max(tpnn$Weight.clust)
tpnn$Weight.normalized <- tpnn$Weight.clust - tpnn$Weight.bp
tpnn$Combined.Weight <- tpnn$Weight.clust + tpnn$Weight.bp
tpnn.no.bp <- tpnn[which(tpnn$Weight.bp==0), c(1,2,5,6)]
names(tpnn.no.bp)[3:4] <- c("Weight", "interaction")
#
# Do pathway.crosstalk.network
pcnn <- pathway.crosstalk.network
dim(pcnn[which(pcnn$interaction=="cluster evidence"),]) #645709 okay, same as total.pathway.net, tpnn
pcnn.clust <- pcnn[which(pcnn$interaction=="cluster evidence"),]
pcnn.clust$Weight <- pcnn.clust$Weight/max(pcnn.clust$Weight)
pcnn.bp <- pcnn[which(pcnn$interaction=="pathway Jaccard similarity"),]
pcnn <- rbind(pcnn.clust, pcnn.bp)
# 
####################################################################################################################
# Question: what is the relationship between GO terms and clusters in bioplanet pathways?
# Cuneyet's file
bp.GO <- readRDS(file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/P2PSim.rds", sep=""))
names(bp.GO)[1:2] <- c("source", "target")
tpn.GO <- merge(tpnn, bp.GO, by=c("source", "target"), all=TRUE)
dev.new()
plot(tpn.GO$Weight.clust ~ tpn.GO$GoBio, pch=19, col=alpha("blue", 0.25))
summary(lm(tpn.GO$Weight.clust ~ tpn.GO$GoBio))
# R-squared:  0.1336
plot(tpn.GO$Weight.clust ~ tpn.GO$GoCel, pch=19, col=alpha("darkgreen", 0.25))
summary(lm(tpn.GO$Weight.clust ~ tpn.GO$GoCel))
# R-squared: 0.1145
plot(tpn.GO$Weight.clust ~ tpn.GO$GoMol, pch=19, col=alpha("red", 0.25))
summary(lm(tpn.GO$Weight.clust ~ tpn.GO$GoMol))
# R-squared: 0.07949
plot(tpn.GO$Weight.bp ~ tpn.GO$GoBio, pch=19, col=alpha("darkblue", 0.25))
summary(lm(tpn.GO$Weight.bp ~ tpn.GO$GoBio))
# R-squared:  0.2208 
plot(tpn.GO$Weight.bp ~ tpn.GO$GoCel, pch=19, col=alpha("turquoise4", 0.25))
summary(lm(tpn.GO$Weight.bp ~ tpn.GO$GoCel))
# R-squared:  0.04385  
plot(tpn.GO$Weight.bp ~ tpn.GO$GoMol, pch=19, col=alpha("magenta", 0.25))
summary(lm(tpn.GO$Weight.bp ~ tpn.GO$GoMol))
# R-squared:  0.1057 
# Contrast GO relationships
plot(tpn.GO$GoCel ~ tpn.GO$GoMol, pch=19, col=alpha("magenta2", 0.25))
summary(lm(tpn.GO$GoCel ~ tpn.GO$GoMol))
# R-squared:0.4674
plot(tpn.GO$GoCel ~ tpn.GO$GoBiol, pch=19, col=alpha("orange3", 0.25))
summary(lm(tpn.GO$GoCel ~ tpn.GO$GoBio))
# R-squared:0.4384 
plot(tpn.GO$GoMol ~ tpn.GO$GoBio, pch=19, col=alpha("green3", 0.25))
summary(lm(tpn.GO$GoMol ~ tpn.GO$GoBio))
# R-squared: 0.6049 
# What pathways are related by GO terms AND clustering? 
test.bio <- tpn.GO[which(tpn.GO$GoBio>0.5 & tpn.GO$Weight.clust>0.5),]
test.cel <- tpn.GO[which(tpn.GO$GoCel>0.5 & tpn.GO$Weight.clust>0.5),]
test.mol <- tpn.GO[which(tpn.GO$GoMol>0.5 & tpn.GO$Weight.clust>0.5),]

save(bioplanet,  bioplanet.list.common, bioplanet.matrix.common, bioplanet.matrix.outersect, bioplanet.matrix.union, bioplanet.adj.matrix, bioplanet.jaccard.matrix, bioplanetgeneedges, bioplanetreledges, bioplanetjaccard.g, bioplanetgenes, bioplanetjaccardedges, get.all.gene.names.from.peps, ambig.gene.clist, ambigs.clist.bioplanet.common, count.ambiguous.gene.weights, count.ambiguous.genes, gene.clist.bioplanet.common.no.ambigs, gene.clist.no.ambigs, get.ambiguous.genes, matrix.common, matrix.union, weighted.matrix.common, matrix.outersect, gene.clist.bioplanet.weighted, gene.clist.weights, pep.clist, gene.clist, calculate.gene.weights.using.pathway, gene.clist.pathway.weights, gene.clist.bioplanet.pathway.weighted, gene.clist.bioplanet.pathway.clustsize.weighted, cluster.pathway.evidence, create.pathway.network, calc.net.density, pathway.net.list, density.list, pathway.net.list2, density.list2, net.atts.df, filter.pathway.edges, pathway.net.list3, big.path.net, pathway.net.list4, thresh.vec, thresh.vec2, zymes.list, tpnn, pcnn, tpn.GO, file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "BioPlanetNetworks.RData", sep=""))


##################################################################################################################
# Focus on enzymes 
# Use Guolin's list of enzymes
enzymes.filename <- paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "Human enzymes_Guolin_06192020.txt", sep="")
# Note: had to delete protein description, columns with links, and other detailed columns for the tab-delimited text file to read. 
zymelist <- read.table(enzymes.filename, header=TRUE, sep = "\t", comment.char = "", na.strings='', stringsAsFactors=FALSE, fill=TRUE)
# Get the primary gene names; separate according to function
names(zymelist)[4] <- "Gene.Name"
zymes.list <- dlply(zymelist, .(Enzyme.Classification.detailed.full), function (x) return (x[,4]))
# 29 enzyme classes
# Check previous classes defined
intersect(ykinases, zymes.list[["Protein-tyrosine kinases"]])
outersect(ykinases, zymes.list[["Protein-tyrosine kinases"]])
# ten differences.
intersect(acetyltransferases, zymes.list[["Acetyltransferase"]])
outersect(acetyltransferases, zymes.list[["Acetyltransferase"]])
# 42 in common; 38 different. E.g. ACAA1 acetyl-CoA acyltransferase 1  is in my older list.
# 
# Can now use this list to look at enrichment of certain enzymes with .common functions above.
#...
zymes.bioplanet.common <- list.common(zymes.list, bioplanet, keeplength = 1)
# This a list of intersecting enzymes and bioplanet pathways 
zymes.gzclusters.common <- list.common(zymes.list, gene.clist, keeplength = 1)
# This a list of intersecting enzymes and clusters from GZ data
# What is the intersection of these two lists?
zymes.gzclusters.bioplanet.common <- list.common(zymes.gzclusters.common, zymes.bioplanet.common, keeplength = 1)
# Write this file for inspection
zymes.crosstalk.filename <- paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "ClusterBioplanetEnzymeIntersect.txt", sep="")
sink(zymes.crosstalk.filename); print(zymes.gzclusters.bioplanet.common); sink()
# What are the ubiquitining enzymes
zymes.list[grep("biquit", names(zymes.list))]
# What enzymes are most represented?
gene.clist["38.104.92"] # FYN YES1 SRC PTPRK PLCG1
# Some more interesting pathways to investigate
bioplanet["CARM1 transcriptional regulation by protein methylation"]
bioplanet["ERBB1 downstream pathway"]
filter.edges.1("Endocytosis", selected.pathway.edges)
intersect(bioplanet[["Endocytosis"]], bioplanet[["Developmental biology"]]) 
endo.paths <- filter.edges.1("Endocytosis", big.path.net)
endo.paths <- endo.paths[endo.paths$Weight > 0.5,]
filter.edges.1("CARM1 transcriptional regulation by protein methylation", big.path.net) 
# Not in selected pathways; max weight = 0.1
erbb1.edges <- filter.edges.1("ERBB1 downstream pathway", big.path.net)
# Not in selected pathway edges; max weight = 1.3
erbb1.edges[erbb1.edges$Weight > 0.5, ]
# IL2 and TGF-beta, FSH, prostaglandin bioshynth; ER protein processing; metabolism; EGF/EGFR siganling pathwyay

# zymes.bioplanet.matrix.common <- matrix.common(bioplanet, bioplanet)

##$####################: 
#Workaround for edge mapping problem in Cytoscape
FixEdgeDprops.RCy32 <- function() {
  allstyles <- getVisualStyleNames()
  setBackgroundColorDefault(style.name = allstyles, "#949494") # grey 58
  setEdgeLineWidthDefault (3. , style.name = allstyles)
  setEdgeColorDefault ( "#FFFFFF", style.name = allstyles)  # white
  setEdgeSelectionColorDefault (col2hex("chartreuse"), style.name = allstyles)
  edgecolors <- col2hex(c("red", "red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2", "darkorange1"))
  edgecolorsplus <- col2hex(c("deeppink", "red", "red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2", "orangered2", "darkorange1"))
  #  red; turquois; green; magenta; blue; violet; green;  bluegreen; black; gray; turquoiseblue; orange 
  edgeTypes <- c("pp", "PHOSPHORYLATION", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of",  "controls-state-change-of", "Physical interactions", "BioPlex", "in-complex-with",  'experiments',  'database',   "Pathway", "Predicted", "Genetic interactions", "correlation", "negative correlation", "positive correlation",  'combined_score', "merged" , "intersect", "peptide", 'homology', "Shared protein domains", "ACETYLATION")
  # 22 edgeTypes            
  myarrows <- c ('Arrow', 'Arrow', 'Arrow','Arrow', 'Arrow', "Arrow", 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', "Arrow")
  setEdgeTargetArrowMapping(style.name = allstyles, 'shared interaction', edgeTypes, myarrows, default.shape='None')  
  # matchArrowColorToEdge('TRUE', style.name = allstyles) # Gives error at this time
  setEdgeColorMapping( 'shared interaction', edgeTypes, edgecolors, 'd', default.color="#FFFFFF", style.name = allstyles )
  # A grey background helps highlight some of the edges
  edgevalues <- getTableColumns('edge')
  if (length(edgevalues[grep("pp", edgevalues$interaction), 1])>0) {
    setEdgeColorBypass(edgevalues[grep("pp", edgevalues$interaction), 1], col2hex("red"))}
  if (length(edgevalues[grep("PHOSPHORYLATION", edgevalues$interaction), 1])>0) {
    setEdgeColorBypass(edgevalues[grep("PHOSPHORYLATION", edgevalues$interaction), 1], col2hex("red"))}
  if (length(edgevalues[grep("phosphorylation", edgevalues$interaction), 1])>0) {
    setEdgeColorBypass(edgevalues[grep("phosphorylation", edgevalues$interaction), 1], col2hex("red"))}
  if (length(edgevalues[grep("ACETYLATION", edgevalues$interaction), 1])>0) {
    setEdgeColorBypass(edgevalues[grep("ACETYLATION", edgevalues$interaction), 1], col2hex("darkorange1"))}
  }
#  
#_________________________________________#_________________________________________
# gzallt.physical.network =	edge file of the physical interaction network: CCCN, CFN, plus negative corr edges
# The simplest way to get the node attributes is to select the genes in the entire network above and the node attribute file.
genenames <- extract.gene.names.RCy3(alkep300)
alkep300.cf <- gz.cf[gz.cf$Gene.Name %in% genenames,]
# Note: this contains PTMs. If you want only the CFN (gene nodes) use:
alkep300.gene.cf <- alkep300.cf[which(alkep300.cf$Node.ID=="gene"),]

# overlap between ldgenes and gzallt genes >>
gzphysicalgenes <- extract.gene.names.RCy3(gzalltgene.physical.cfn.merged)
ldphysicalgenes <- extract.gene.names.RCy3(ldgene.physical.cfn.merged)
setdiff(gzphysicalgenes, ldphysicalgenes)
setdiff(ldphysicalgenes, gzphysicalgenes)
outersect(ldphysicalgenes, gzphysicalgenes)
setdiff(V(gzalltgene.physical.cfn.merged.g),V(ldgene.physical.cfn.merged.g))
setdiff(V(ldgene.physical.cfn.merged.g),V(gzalltgene.physical.cfn.merged.g))
head(V(ldgene.physical.cfn.merged.g))
head(V(gzalltgene.physical.cfn.merged.g))
length(ldphysicalgenes)
length(gzphysicalgenes)
class(V(gzalltgene.physical.cfn.merged.g))
gzigraphgenes <- unlist(V(gzalltgene.physical.cfn.merged.g))
ldigraphgenes <- unlist(V(ldgene.physical.cfn.merged.g))
as.character(ldigraphgenes)
as.vector(ldigraphgenes)
class(ldigraphgenes)
# ?V: A vertex sequence is just what the name says it is: a sequence of vertices.Return: A vertex sequence containing all vertices, in the order of their numeric vertex ids.
difference(gzigraphgenes, ldigraphgenes)
difference(ldigraphgenes, gzigraphgenes)
gzva <- unlist(vertex_attr(gzalltgene.physical.cfn.merged.g))
ldva <- unlist(vertex_attr(ldgene.physical.cfn.merged.g))
setdiff(gzva, ldva)
setdiff(ldva, gzva)
#*#*#*#*
#* Hypothesis: clustering mulitiple edges will reveal patterns in interactions
# * Plot multiple edge weights in 2D or 3D and use minimum spanning tree, single linkage
# * For more than 3D, use t-SNE to find patterns
#   * LF Data: Gene CCCN, CFN; PTM Gene CCCN, CFN, bioplanet jaccard similarity
#   * Define bioplanet jaccard similarity between genes as that for the pathway interactions
#
# Use approach for making CCCN from clusterlist
#   make.clusterlist is function(tsnedata, toolong, tbl.sc)
#   Try with commandes that are built into that function
tpn.weights <- total.pathway.net[, c("Weight.bp",  "Weight.clust")]
# plot a subset of the 645709 edges for testing
tpnw.sub <- tpn.weights[seq(1, 645000, 100), ]
plot(tpnw.sub)
plot(log2(tpnw.sub))
# Convert to log2
tpnw.sub.log2 <- log2(tpnw.sub)
tpnw.log2 <- log2(tpn.weights)
plot(tpnw.log2)
toolong=0.1
tpn.span2 <- spantree(dist(tpnw.sub), toolong=toolong)
tpn.disc2 <-  distconnected(dist(tpnw.sub), toolong = toolong, trace = TRUE)  # test
cat ("threshold dissimilarity", toolong, "\n", max(tsnedata.disc2), " groups","\n")
ordiplot(tpnw.sub)
plot(tpnw.sub)
plot(tpnw.log2)
# Conclusion: there is not enough strucure in the data -- it's a cloud, very dense!
#lines(tsne.span2, tsnedata)
ordihull(tpnw.sub.log2, tpn.disc2, col="red", lwd=2)	
# Find groups
tsnedata.span2.df <- data.frame(rownames(tbl.sc))
names(tsnedata.span2.df) <- "Gene.Name"
tsnedata.span2.df$group <- tsnedata.disc2
tsnedata.span2.list <- dlply(tsnedata.span2.df, .(group))  # GROUP LIST  !
return(tsnedata.span2.list)	



