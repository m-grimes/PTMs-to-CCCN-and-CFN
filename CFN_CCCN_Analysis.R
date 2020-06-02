###############
#Functions for shortest paths analysis
##############

#Select sites with absolute value greater than a specified minimum
select.sites.w.min.value <- function(df, exp, min_value = 2.25) {
  data <- df[, exp]
  df <- df[ which(abs(data) > log2(min_value)), ]
  return(df)
    
}

#Select sites with at least the desired number of observations (not NA); default minimum observations is 2)
select.sites.w.enough.obs <- function(df, num_rep, min_obs = 2) {
  zero.to.NA(df)
  df$sumNA <- apply(is.na(df), 1, sum)
  allowed.na <- num_rep - min_obs
  df <- df[ which(df$sumNA <= allowed.na), ]
  return(df)
}

#Select sites with absolute value of mean at least the specific value.
#This function assumes that the columns with the ratio are log2 transformed and consecutive
#The first_rep parameter is the number of the first column with data

select.sites.w.min.mean <- function(df, first_rep = 3, num_rep, min_mean = 2.25) {
  last_rep = first_rep + num_rep - 1
  df$mean <- rowMeans(df[first_rep:last_rep], na.rm = TRUE)
  df <- df[ which(abs(df$mean) > log2(min_mean)), ]
  return(df)
  
}

#Select sites with coefficient of variation (log-transformed) less than or equal to the specified value.
#This function assumes that the columns with the ratio are log2 transformed and consecutive
#The first_rep parameter is the number of the first column with data
select.sites.w.max.cv <- function(df, first_rep = 3, num_rep, max_cv = 1) {
  last_rep = first_rep + num_rep - 1
  df$sd_base2 <- apply(df[first_rep:last_rep],1, sd, na.rm = TRUE)
  df$sd_ln <- df$sd_base2 * log(2)
  df$cv <- sqrt(exp(df$sd_ln^2) -1)
  df <- subset(df, cv <= 1)
  return(df)
}

#Create the shortest paths sub-network of the CFN that connects a list of sites to a specified target.
#This functions assumes a dataframe with a row for each selected site and a column called Gene.Name with the gene symbol for each site
shortest.paths.to.target <- function(df, target) {
  df_changed_genes <- unique(df$Gene.Name)
  df_changed_genes2 <- df_changed_genes[df_changed_genes %in% gzallt.gene.key$Gene.Name]
  df_paths <- composite.shortest.paths(genes1=c(target), df_changed_genes2, exclude="", network=gzalltgene.physical.cfn.merged)
  return(df_paths)
}

#Takes any number of dataframes as arguments and returns the union of the dataframes; 
#Uses the bind_rows function from the dplyr package
my.union <- function(...) {
    dfs <- list(...)
    dfs_union <- bind_rows(dfs)
    dfs_union <- dfs_union[!duplicated(dfs_union),]
    return(dfs_union)
}

#Takes any number of dataframes as arguments and returns the intersection
my.intersect <- function(...) {
  dfs <- list(...)
  dfs_intersect <- Reduce(intersect, dfs)
  return(dfs_intersect)
}

