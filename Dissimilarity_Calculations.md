Dissimilarity Representation Calculations
================
Mark Grimes, Ekaterina Smirnova
1/30/2020

## R Markdown

``` r
require(rgl)
```

    ## Loading required package: rgl

``` r
require(Rtsne)
```

    ## Loading required package: Rtsne

``` r
require(vegan)
```

    ## Loading required package: vegan

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-6

``` r
require(pryr)
```

    ## Loading required package: pryr

    ## Registered S3 method overwritten by 'pryr':
    ##   method      from
    ##   print.bytes Rcpp

``` r
require(plyr)
```

    ## Loading required package: plyr

``` r
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::arrange()   masks plyr::arrange()
    ## x purrr::compact()   masks plyr::compact()
    ## x purrr::compose()   masks pryr::compose()
    ## x dplyr::count()     masks plyr::count()
    ## x dplyr::failwith()  masks plyr::failwith()
    ## x dplyr::filter()    masks stats::filter()
    ## x dplyr::id()        masks plyr::id()
    ## x dplyr::lag()       masks stats::lag()
    ## x dplyr::mutate()    masks plyr::mutate()
    ## x purrr::partial()   masks pryr::partial()
    ## x dplyr::rename()    masks plyr::rename()
    ## x dplyr::summarise() masks plyr::summarise()
    ## x dplyr::summarize() masks plyr::summarize()

``` r
# Focus on intersect of all clusters from Euclid, Spearman, and SED
list.common <- function (list1, list2, keeplength=3) {
  parse <- lapply(list1, function (y) sapply(list2,  function(x) intersect(x, y)))
  dims <- lapply(parse, function (x) sapply(x, length))
  keep <- which(sapply(dims, sum) > keeplength)
  pare <- parse[keep]
  prune <- lapply(pare, function (y) return (y[which(sapply(y, function (x) which(length(x) > keeplength )) > 0)]))
  newlist <- unlist(prune, recursive=FALSE)
  return(newlist)
}

# Define clusters from t-SNE embeddings
make.clusterlist <- function(tsnedata, toolong, tbl.sc, title = "") {
  tsne.span2 <- spantree(dist(tsnedata), toolong=toolong)
  tsnedata.disc2 <-  distconnected(dist(tsnedata), toolong = toolong, trace = FALSE)  # test
  cat ("threshold dissimilarity", toolong, "\n", max(tsnedata.disc2), " groups","\n")
  p1.pryr %<a-% {
    par(mar=c(1,1,1,1))
    ordiplot(tsnedata, main = title)
    ordihull(tsnedata, tsnedata.disc2, col="red", lwd=2)
  }
  # Find groups
  tsnedata.span2.df <- data.frame(rownames(tbl.sc))
  names(tsnedata.span2.df) <- "Peptide.Name"
  tsnedata.span2.df$group <- tsnedata.disc2
  tsnedata.span2.list <- dlply(tsnedata.span2.df, .(group))  # GROUP LIST  !
  return(list(tsnedata.list = tsnedata.span2.list, p1.pryr))    
}

make.adj.mat <- function(list.element)      {
  list.el.mat <- matrix(1, nrow=length(list.element), ncol=length(list.element))
  rownames(list.el.mat) <- list.element
  colnames(list.el.mat) <- list.element
  return(list.el.mat)
}
extract.genes.from.clist <- function (clusterlist.element) {
  element <- clusterlist.element[1]
  genes <- unique(sapply(as.character(element$Peptide.Name),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
  return(genes)
}
extract.peps.from.clist <- function (clusterlist.element) {
  element <- clusterlist.element[1]
  return(as.character(element$Peptide.Name))
}
zero.to.NA <- function(df) {
  zer0 <- which(df==0, arr.ind = TRUE)
  cfNA <- as.matrix(df)
  cfNA[zer0] <- NA
  cfNA <- data.frame(cfNA)
  return(cfNA)
}
```

## Dissimilarity and distance calculations

This is Katia’s version:

Katia’s experimental functions to generate embeddings

``` r
#Utilize calc.dist function to perform tsne clustering
calc.tsne <- function(ptm.data, plots_path_dir, dims = 3, perplexity = 15, theta = 0.25, 
                      check_duplicates = FALSE, pca=FALSE, toolong = 3.5, radius=0.76,
                      nbreaks.sizes = 100, nbreaks.common.sizes = 200){
  
  #calculate each of the 3 PTM distance matrices using calc.dist function
  #euclidean distance
  eucl.dist <- calc.dist(ptm.data = ptm.data, dist.type = "euclidean")
  #spearman dissimilarity
  spearman.dist <- calc.dist(ptm.data = ptm.data, dist.type = "spearman.dissim")
  #a combination of spearman and euclidean distances
  sed.dist <- calc.dist(ptm.data = ptm.data, dist.type = "sed", 
                       cor.mat = list(spearman.mat = spearman.dist$cor.mat$spearman.mat, eucl.mat = eucl.dist$cor.mat$eucl.mat))
  #tsne clusters using euclidean distance
  eucl.tsne.list <- Rtsne(eucl.dist$dist.mat, dims = dims, perplexity = perplexity, theta = theta, 
                          check_duplicates = check_duplicates, pca=pca)
  eucl.tsne <- eucl.tsne.list$Y
  #tsne clusters using spearman dissimilarity
  spearman.tsne.list <- Rtsne(spearman.dist$dist.mat, dims = dims, perplexity = perplexity, theta = theta, 
                          check_duplicates = check_duplicates, pca=pca)
  spearman.tsne <- spearman.tsne.list$Y
  #tsne clusters using sed dissimilarity
  sed.tsne.list <- Rtsne(sed.dist$dist.mat, dims = dims, perplexity = perplexity, theta = theta, 
                              check_duplicates = check_duplicates, pca=pca)
  sed.tsne <- sed.tsne.list$Y
  
  # Use the above embeddings (using gz data by itself) to define clusters
  # Note: in practice the toolong value is adjusted to so that clusters are defined in a way that makese sense empirically. Formally, values around 3.5 are considered a good boundary for t-SNE embeddings.
  pdf(file.path(plots_path_dir, "Hull.pdf"))
  eucl.cluster.list <- make.clusterlist(eucl.tsne, toolong = toolong, 
                        tbl.sc = ptm.data, title = "Euclidean distance clusters")$tsnedata.list
  eucl.cluster.sizes <- sapply(eucl.cluster.list, function(x) dim(x)[1])
  spearman.cluster.list <- make.clusterlist(spearman.tsne, toolong = toolong, 
                              tbl.sc = ptm.data, title = "Spearman distance clusters")$tsnedata.list 
  
  spearman.cluster.sizes <- sapply(spearman.cluster.list, function(x) dim(x)[1])
  sed.cluster.list <- make.clusterlist(sed.tsne, toolong=toolong, 
                      tbl.sc = ptm.data, title = "Spearman and Euclidean distance clusters")$tsnedata.list
  sed.cluster.sizes <- sapply(sed.cluster.list, function(x) dim(x)[1])
  dev.off()
  
  #return cluster characteristics
  # Three dimensional plots
  plot3d(eucl.tsne, type="s", radius=radius, 
         col="forestgreen", xlab= "Dim 1", ylab = "Dim 2", zlab = "Dim 3")
  writeWebGL(dir = plots_path_dir, filename = file.path(plots_path_dir, "Eucl.html"),
             width = 500, reuse = TRUE)
  
  #rgl.snapshot(file.path(plots_path_dir, "Eucl.png"))
  plot3d(spearman.tsne, type="s", 
         radius=radius, col="purple",  xlab= "Dim 1", ylab = "Dim 2", zlab = "Dim 3")
  writeWebGL(dir = plots_path_dir, filename = file.path(plots_path_dir, "Spearman.html"),
             width = 500, reuse = TRUE)
  #rgl.snapshot(file.path(plots_path_dir, "Spearman.png"))
  
  plot3d(sed.tsne, type="s", 
         radius=radius, col="red",  xlab= "Dim 1", ylab = "Dim 2", zlab = "Dim 3") 
  writeWebGL(dir = plots_path_dir, filename = file.path(plots_path_dir, "Sed.html"),
             width = 500, reuse = TRUE)
  #rgl.snapshot(file.path(plots_path_dir, "Sed.png"))
  
  
  pdf(file.path(plots_path_dir, "NClusters.pdf"))
  #histogram of cluster sizes
  par(mfrow = c(3,1))
  hist(eucl.cluster.sizes, breaks=nbreaks.sizes, col="green",
       main = "Euclidean tsne cluster sizes", xlab = "cluster size") 
  # note two or three large clusters that should be broken up
  hist(spearman.cluster.sizes, breaks=nbreaks.sizes, col="green",
       main = "Spearman tsne cluster sizes", xlab = "cluster size") 
  # note 2-4 large clusters that should be broken up
  hist(sed.cluster.sizes, breaks=nbreaks.sizes, col="green",
       main = "Sed tsne cluster sizes", xlab = "cluster size") 
  dev.off()
  
  #the next step is to examine cluster histograms of each distance type 
  #and identify whether they should be broken further
  # note two large clusters that should be broken up
  #  Concatenate group list files
  #     make group names unique
  #
  eucl.cluster <- ldply(eucl.cluster.list)[,2:3]
  spearman.cluster <- ldply(spearman.cluster.list)[,2:3]
  sed.cluster <- ldply(sed.cluster.list)[,2:3]  # Further partition large groups
  eucl.cluster$group <- paste(noquote(eucl.cluster$group), noquote("e"), sep="", collapse=NULL)
  spearman.cluster$group <- paste(noquote(spearman.cluster$group), noquote("s"), sep="", collapse=NULL)
  sed.cluster$group <- paste(noquote(sed.cluster$group), noquote("sed"), sep="", collapse=NULL)
  cluster.df <- rbind(eucl.cluster, spearman.cluster, sed.cluster)
  # 9396
 

  eucl.genes <- lapply(eucl.cluster.list, extract.genes.from.clist)
  spearman.genes <- lapply(spearman.cluster.list, extract.genes.from.clist)
  sed.genes <- lapply(sed.cluster.list, extract.genes.from.clist)
  #
  eucl.peps <- lapply(eucl.cluster.list, extract.peps.from.clist)
  spearman.peps <- lapply(spearman.cluster.list, extract.peps.from.clist)
  sed.peps <- lapply(sed.cluster.list, extract.peps.from.clist)
  

  eu.sp.gz <- list.common(eucl.peps, spearman.peps, keeplength=2)
  eu.sp.gz.sizes <- sapply(eu.sp.gz, length)
  common.peps.cluster <- list.common(eu.sp.gz, sed.peps, keeplength=2)#defines overlap and require its min 3 by default
  common.peps.cluster.sizes <- sapply(common.peps.cluster, length) 
  
  #histogram of common cluster sizes
  pdf(file.path(plots_path_dir, "CommonPepsSizes.pdf"))
  par(mfrow = c(1,1))
  hist(common.peps.cluster.sizes, breaks=nbreaks.common.sizes, col="gold",
       main = "Common tsne cluster sizes", xlab = "cluster size")
  dev.off()
  
  #list of the 3 dissmimilarity matrices
  dist.mat.list <- list(eucl.dist = eucl.dist, spearman.dist = spearman.dist, sed.dist = sed.dist)
  #list of the 3 tsne cluster results
  cluster.list <- list(eucl.cluster.list = eucl.cluster.list, spearman.cluster.list = spearman.cluster.list,
       sed.cluster.list = sed.cluster.list)
  
  #list of cluster sizes for plotting each distance histogram
  cluster.sizes <- list(eucl.cluster.sizes = eucl.cluster.sizes, 
                        spearman.cluster.sizes =spearman.cluster.sizes,
                        sed.cluster.sizes =sed.cluster.sizes)
  
  #list of the 3 tsne gene results that belong to each cluster
  genes.clusters <- list(eucl.genes = eucl.genes, spearman.genes = spearman.genes, 
                         sed.genes = sed.genes) 
  
  #list of the 3 tsne peptites results that belong to each cluster
  peps.clusters <- list(eucl.peps = eucl.peps, spearman.peps = spearman.peps, 
                         sed.peps = sed.peps)
  
  #############################################
  #construct adjecency matrix based on co-clustered peptites
  ##############################################

  # remove self loops
  diag(spearman.dist$dist.mat) <- NA
  #create an edge (put 1 in adjacency matrix) if the two peptites are in the same cluster
  #based on the peptites that end up in the same cluster using tsne based on 3 distance metrics
  gz.adj <- rbind.fill.matrix(llply(common.peps.cluster, make.adj.mat))
  rownames(gz.adj) <- colnames(gz.adj)
  #dim: 2764 2764
  gz.adj.o <- gz.adj[order(rownames(gz.adj)), order(colnames(gz.adj))]
  #select correlation values for the  peptites 
  #that appear in all clusters using tsne based on 3 distance metrics 
  gz.cccn <- spearman.dist$dist.mat[intersect(rownames(gz.adj.o), rownames(spearman.dist$dist.mat)), intersect(colnames(gz.adj.o), colnames(spearman.dist$dist.mat))]
  #gz.cccn.1  and gz.cccn are identical, why two matrices?
  gz.NA <- which(is.na(gz.adj.o), arr.ind = TRUE)
  #fill in the na in cluster correlated network if it is an NA in the adjecency matrix
  gz.cccn <- replace (gz.cccn, gz.NA, NA)
  # remove self loops
  diag(gz.cccn) <- NA
  
  return(list(dist.mat.list = dist.mat.list,
              cluster.list = cluster.list, cluster.sizes = cluster.sizes,
              cluster.df = cluster.df, common.peps.cluster = common.peps.cluster,
              genes.clusters = genes.clusters, peps.clusters = peps.clusters,
              cluster.network = gz.cccn))
}

gene.adj.mat <- function(cluster.network,  cor.threshold =0.5){
  gz.cccn = cluster.network#change naming conventions later
  ## ****
  # Option: Limit to a particlar correlation value like this.
  gz.cccn.halflim <- replace (gz.cccn, abs(gz.cccn) < cor.threshold, NA)
  ##################################################################
  # Construct CCCN for proteins (genes) based on correlations of mod sites: 
  # Make Gene CCCN
  # Note: Did this on server for large data set, much faster at ddply step
  gz.gene.cccn <- data.frame(gz.cccn, row.names = rownames(gz.cccn), 
                             check.rows=TRUE, check.names=FALSE, 
                             fix.empty.names = FALSE)
  gz.gene.cccn$Gene.Name <- sapply(rownames(gz.gene.cccn), function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])
  #this is smaller than in the KarenGuolinData3.R file
  length(unique(gz.gene.cccn$Gene.Name))   # 1346 for first test
  # Use only upper triangle so correlations are not duplicated during the next step
  gz.gene.cccn[lower.tri(gz.gene.cccn)] <- NA   
  # ________________________________
  # Sum correlations in one dimension, then the other dimension
  #what is this doing??? why you need it?
  #2764 is the number of peptites we take correlations on,
  #the correlations on those peptites are groupped on the gene level 
  #in the columns of gz.gene.cccn2  -- don't understand why sum correlations? 
  #what if correlations sum > 1?
  gz.gene.cccn2 <- ddply(gz.gene.cccn, .(Gene.Name), numcolwise(function(x) sum(abs(x), na.rm=T)), .progress = "tk")
  dim(gz.gene.cccn2)    #   1429 2765 -- mine are smaller (2764 1347), gene names changed???
  rownames(gz.gene.cccn2) <- gz.gene.cccn2$Gene.Name
  gz.gene.cccn2 <- gz.gene.cccn2[, 2:ncol(gz.gene.cccn2)]
  gz.gene.cccn2 <- data.frame(t(gz.gene.cccn2))
  gz.gene.cccn2$Gene <- sapply(rownames(gz.gene.cccn2), function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])
  # other dimension -- again sum correlations across peptites (in rows) 
  #to have them on th egene level
  gz.gene.cccn3 <- ddply(gz.gene.cccn2, .(Gene), numcolwise(function(x) sum(x, na.rm=T)), .progress = "tk")
  dim(gz.gene.cccn3)  #  1429 1430  (mine are 1346 1347 instead)
  
  ####################################
  #fix fix part to make it more reliable 
  #in case there are more gene names with problematic namein conventions
  #in the gene names to avoid checking it for every data set
  #the gene names have '-', which is later replaced by a '.' if the data is converted into 
  #data frame
  #####################################
  # R likes to put dots in column names, which is a problem for ambiguous gene names and gene names with hyphens 
  # E.g., Note punctuation missing for "HLA"  "NKX2" "NME1 in column names
  gz.gene.cccn3$Gene[grep("HLA", gz.gene.cccn3$Gene)]
  names(gz.gene.cccn3)[grep("HLA", names(gz.gene.cccn3))]
  fix.baddot <- function(cell) {  
    baddot <-  c("HLA.A", "HLA.B","HLA.C","HLA.F", "HLA.H", "NKX2.1","NME1.NME2")
    goodhyphen <-  c("HLA-A","HLA-B","HLA-C","HLA-F","HLA-H","NKX2-1","NME1-NME2")
    cellv <- unlist(strsplit(as.character(cell), "; "))
    if (any(baddot %in% cellv)) {
      cellv.new <- gsub(baddot[baddot %in% cellv], goodhyphen[baddot %in% cellv], cellv)    
      return (paste(cellv.new, collapse="; "))
    } else return(cell)     }
  #
  fixednames <- sapply(names(gz.gene.cccn3), fix.baddot)    
  fixednames[grep("HLA", fixednames)]
  # Also if necessary
  fixednames <- sapply(fixednames, function (x) gsub("\\.", ";", x))
  
  # Okay this is a pain, so just work around the problem (once satisfied that the gene names actually match). 
  names(gz.gene.cccn3)[2:ncol(gz.gene.cccn3)] <- gz.gene.cccn3$Gene
  rownames(gz.gene.cccn3) <- gz.gene.cccn3$Gene
  identical (rownames(gz.gene.cccn3), names(gz.gene.cccn3[,2:ncol(gz.gene.cccn3)])) # TRUE
  
  # Make a adjacency matrix with NAs to index edges only for genes that co-cluster
  gz.gene.cccn0 <- gz.gene.cccn3[,2:ncol(gz.gene.cccn3)]

  gz.gene.cccn.na <- zero.to.NA(gz.gene.cccn0)
  # Make igraph object
  gz.gene.cccn.g <- graph.adjacency(as.matrix(gz.gene.cccn0), mode="lower", diag=FALSE, weighted="Weight")
  # Here is how to make an edge list file
  gzgenecccn.edges <- data.frame(as_edgelist(gz.gene.cccn.g))
  names(gzgenecccn.edges) <- c("source", "target")
  gzgenecccn.edges$Weight <- edge_attr(gz.gene.cccn.g)[[1]]
  gzgenecccn.edges$interaction <- "correlation" 
  gzgenecccn.edges $interaction[gzgenecccn.edges$Weight<=-0.5] <- "negative correlation"
  gzgenecccn.edges $interaction[gzgenecccn.edges$Weight>=0.5] <- "positive correlation"
  # The following is useful for driving cytoscape layouts to push negative correlation edges away (Cytoscape doesn't recognize negative values for edge weight).
  negs <- gzgenecccn.edges[gzgenecccn.edges$Weight<0,'Weight']
  frnegs <- exp (negs*20)
  
  gzgenecccn.edges$Alt.Weight <- gzgenecccn.edges $Weight
  gzgenecccn.edges[gzgenecccn.edges$Weight<0,'Alt.Weight'] <- frnegs    
  gzgenecccn.edges[grep("SRC", gzgenecccn.edges$source),]
  
  return(gzgenecccn.edges = gzgenecccn.edges)
  
}#end of function that creates adjacency matrix
```
