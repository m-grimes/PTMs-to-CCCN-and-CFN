require(rgl)
require(Rtsne)
require(vegan)
require(plyr)
#function to 
#(1) calculate correlation matrix on PTMs
#(2) use  correlation matrix to construct tsne embeddings

#input: ptm.data = gzdata - a data frame with PTMs in rows and samples in columns
#dist.type: c("euclidean", "spearman.dissim", "sed")
#cor.mat: - if method = "dissimilarity" user can provide spearman correlation matrix 
#to avoid computing spearman correlation matrix twice

calc.dist <- function(ptm.data, dist.type, cor.mat = list(spearman.mat = NULL, eucl.mat = NULL)){
  
  # set NA to two orders of magnitude higher than max distance
  if(dist.type == "euclidean"){
    eucl.mat = as.matrix (dist (ptm.data), method = "euclidean")
    dist.mat = eucl.mat
    dist.mat[is.na(dist.mat)] <- 100*max(dist.mat, na.rm=T)
    cor.mat$eucl.mat = eucl.mat}
  else if(dist.type == "spearman.dissim"){
    spearman.mat <- cor(t(ptm.data), use = "pairwise.complete.obs", method = "spearman")
    dist.mat <- 1 - abs(spearman.mat)
    # set NA to two orders of magnitude higher than max distance
    dist.mat[is.na(dist.mat)] <- 100*max(dist.mat, na.rm=T)
    cor.mat$spearman.mat = spearman.mat
  }
  else if(dist.type == "sed"){
    #combine Euclid and Spearman w/o taking absolute value.	
    spearman.mat <- cor.mat$spearman.mat
    eucl.mat <- cor.mat$eucl.mat
    if(is.null(spearman.mat)){spearman.mat <- cor(t(ptm.data), use = "pairwise.complete.obs", method = "spearman")}
    diss.noabs <- 1 - spearman.mat  
    diss.noabs[is.na(diss.noabs)] <- 50*max(diss.noabs, na.rm=T) 
    
    if(is.null(eucl.mat)){eucl.mat = as.matrix (dist (ptm.data), method = "euclidean") }
    
    eucl.mat[is.na(eucl.mat)] <- 100*max(eucl.mat, na.rm=T)
    eucl.mat.1 <- 100*eucl.mat/max(eucl.mat, na.rm=T)  # now max=100
    # SED: combine Euclid and Spearman w/o taking absolute value.
    dist.mat <- (eucl.mat.1 + diss.noabs)/2
  }
  
  return(list(dist.mat = as.matrix(dist.mat), cor.mat =cor.mat))
}

#Utilize calc.dist function to perform tsne clustering
calc.tsne <- function(ptm.data, plots_path_dir, dims = 3, perplexity = 15, theta = 0.25, 
                      check_duplicates = FALSE, pca=FALSE, toolong = 3.5, radius=0.76){
  
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
  
  
  # Define clusters from t-SNE embeddings
  make.clusterlist <- function(tsnedata, toolong, tbl.sc, title = "")	{
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
  #browser()
  #return cluster characteristics
  # Three dimensional plots
  #writeWebGL(dir = plots_path_dir, filename = file.path(plots_path_dir, "Eucl.html"),
  #           width = 500, reuse = TRUE)
  plot3d(eucl.tsne, type="s", radius=radius, 
         col="forestgreen", xlab= "Dim 1", ylab = "Dim 2", zlab = "Dim 3")
  rgl.snapshot(file.path(plots_path_dir, "Eucl.png"))
  #writeWebGL(dir = plots_path_dir, filename = file.path(plots_path_dir, "Spearman.html"),
  #           width = 500, reuse = TRUE)
  plot3d(spearman.tsne, type="s", 
         radius=radius, col="purple",  xlab= "Dim 1", ylab = "Dim 2", zlab = "Dim 3")
  rgl.snapshot(file.path(plots_path_dir, "Spearman.png"))
  #writeWebGL(dir = plots_path_dir, filename = file.path(plots_path_dir, "Sed.html"),
  #           width = 500, reuse = TRUE)
  plot3d(sed.tsne, type="s", 
         radius=radius, col="red",  xlab= "Dim 1", ylab = "Dim 2", zlab = "Dim 3") 
  rgl.snapshot(file.path(plots_path_dir, "Sed.png"))
  
  
  pdf(file.path(plots_path_dir, "Hist.pdf"))
  #histogram of cluster sizes
  hist(eucl.cluster.sizes, breaks=100, col="green") 
  # note two or three large clusters that should be broken up
  hist(spearman.cluster.sizes, breaks=100, col="green") 
  # note 2-4 large clusters that should be broken up
  hist(sed.cluster.sizes, breaks=100, col="green") 
  dev.off()
  
  return(list(eucl.cluster.list = eucl.cluster.list, spearman.cluster.list = spearman.cluster.list,
              sed.cluster.list = sed.cluster.list))
}

