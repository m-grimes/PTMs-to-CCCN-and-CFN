# BioPlanet Pathway Interaction Networks
# July 9, 2021
# Mark Grimes

# =============================================================================
# Construct a network of bioplanet pathway relationships based on lung cancer PTM data
# Summary of steps created and tested in "BioPlanet_Networks.R" ______________________________________________________________________________
#==============================================     
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
# How many pathways per gene?
calculate.pathways.per.gene <- function(gene, pathwaylist) {
  sum.paths <- sum(unlist(sapply(pathwaylist, function (x) which(gene %in% x))))
  return(sum.paths)
}
bioplanetgenes$no.pathways <- sapply(bioplanetgenes$Gene.Name, calculate.pathways.per.gene, pathwaylist=bioplanet)
# This still took over an hour!
hist(as.numeric(bioplanetgenes$no.pathways), breaks=100, col="aquamarine")
head(bioplanetgenes[order(bioplanetgenes$number.of.pathways, decreasing=TRUE),], 30)
# Behold the usual suspects. 
#
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
#
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
##
#----------------------------------------------------------------------

# How does this filter pathway-pathway relationships?
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

# Call edge weight the sum of evidence for interaction
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
# Note: set threshold to zero because many small interactions add up.
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
#
calc.net.density <- function(edgefile) {
  no.nodes <- length(unique(c(edgefile$source, edgefile$target)))
  no.edges <- dim(edgefile)[2]
  no.possible.edges <- 0.5*no.nodes*(no.nodes-1)
  net.density <- no.edges/no.possible.edges
  return(net.density)
}
calc.net.density(total.net) # 3.615564e-06

plot(total.pathway.net$Combined.Weight~total.pathway.net$Weight.clust)
# 456985/645709 edges
head(total.pathway.net.no.bp[1:25, 1:2], 25)
# Metabolism; RNA binding proteins, messengry RNA processing; splicing; feature prominently 
#NOTE: big.path.net filters individual edges to be >0.01; weights are larger in total.pathway.net
total.pathway.net <- total.pathway.net[order(total.pathway.net$Combined.Weight,decreasing=TRUE),]
dev.new()
plot(total.pathway.net$Weight.bp~total.pathway.net$Weight.clust, log="xy", pch=20, col=alpha(col2hex("seagreen4"), 0.1), cex=1.8)
# A cloud
plot(total.pathway.net$Combined.Weight~total.pathway.net$Weight.clust, log="xy", pch=20, col=alpha(col2hex("magenta"), 0.1), cex=1.8)
points(big.path.net$Weight~big.path.net$Weight.clust,  pch=19, col=alpha(col2hex("blue"), 0.1), cex=1.8)
points(total.pathway.net$Weight.normalized~total.pathway.net$Weight.clust,  pch=20, col=alpha(col2hex("yellow"), 0.1), cex=1.8)
###############################################################################
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
# SAVE
bioplanetjaccard.g <- bioplanet.g

save(bioplanet,  bioplanet.list.common, bioplanet.matrix.common, bioplanet.matrix.outersect, bioplanet.matrix.union, bioplanet.adj.matrix, bioplanet.jaccard.matrix, bioplanetgeneedges, bioplanetreledges, bioplanetjaccard.g, bioplanetgenes, bioplanetjaccardedges, get.all.gene.names.from.peps, ambig.gene.clist, ambigs.clist.bioplanet.common, count.ambiguous.gene.weights, count.ambiguous.genes, gene.clist.bioplanet.common.no.ambigs, gene.clist.no.ambigs, get.ambiguous.genes, matrix.common, matrix.union, weighted.matrix.common, matrix.outersect, gene.clist.bioplanet.weighted, gene.clist.weights, pep.clist, gene.clist, calculate.gene.weights.using.pathway, gene.clist.pathway.weights, gene.clist.bioplanet.pathway.weighted, gene.clist.bioplanet.pathway.clustsize.weighted, cluster.pathway.evidence, create.pathway.network, calc.net.density, pathway.net.list, density.list, pathway.net.list2, density.list2, net.atts.df, filter.pathway.edges, pathway.net.list3, big.path.net, pathway.net.list4, thresh.vec, thresh.vec2, zymes.list, create.pathway.network, total.pathway.net, total.pathway.net.no.bp, pathway.crosstalk.network, file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "BioPlanetNetworks.RData", sep=""))
