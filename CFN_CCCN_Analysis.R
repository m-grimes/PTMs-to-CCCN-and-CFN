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
#This functions assumes a dataframe with a row for each selected site and a column called Gene.
#Name with the gene symbol for each site
shortest.paths.to.target <- function(df, target,excluded_gene= "") {
  df_changed_genes <- unique(df$Gene.Name) #unique genes in the dataframe
  df_changed_genes2 <- df_changed_genes[df_changed_genes %in% gzallt.gene.key$Gene.Name] #unique genes in the experimental dataset
  df_paths <- composite.shortest.paths(genes1=c(target), 
                                       df_changed_genes2, exclude=excluded_gene, 
                                       network=gzalltgene.physical.cfn.merged)
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

###############
#Functions for Cluster Enrichment Analysis
##############

#This function takes a dataframe showing which sites are in which clusters, a dataframe of significantly changed sites, 
#and a dataframe with a count of how many sites are in (and not in) each cluster. It maps the sites to clusters, 
#tallies the number of sig sites in each cluster (and not in each cluster), and performs the fisher test. 
#Returns a dataframe with the values used for the fisher 2x2 contingency table, the raw p-value, and the Benjamini 
#corrected p-value.
cluster_enrichment <- function(clusters, sig_sites, all_sites_tally) {
  #Change name of id column of sig_sites table to "Site" so it can be merged with the clusters table
  colnames(sig_sites)[1] = "Site" 
  
  #To map the sig sites to clusters do a merge of the sig_sites to the clusters table. 
  #Note that some of the sig sites do not map to a cluster, so the number of rows in the merged table is less than the number
  #in the sig sites table
  mapped_sites <- merge(clusters, sig_sites, by = "Site")
  
  #Calculate the number of selected sites in each cluster
  selected_tally <- plyr::count(mapped_sites, var = "Cluster")
  colnames(selected_tally)[2] <- "SelectedNumInCluster"
  
  #Calculate the number of selected sites not in each cluster. 
  num_selected_sites <- nrow(mapped_sites)
  selected_tally$SelectedNumOutCluster <- num_selected_sites - selected_tally$SelectedNumInCluster
  
  #Merge the data for all clusters and the data for the selected clusters to make a table with all data needed for the
  #Fisher test. Set the row names to be the cluster ID's. Remove rows where the number of selected sites in a cluster
  #is less than 3 because this is unlikely to give a significant result and will just drive up the number of tests
  fisher_table <- merge(selected_tally, all_tally, by = "Cluster")
  rownames(fisher_table) <- fisher_table$Cluster
  fisher_table$Cluster <- NULL
  fisher_table <- fisher_table[ which(fisher_table$SelectedNumInCluster >= 3),]
  
  #Perform the Fisher Test and the Benjamini-Hochberg multiple testing correction
  fisher_table$p_value <- apply(fisher_table,1, function(x) fisher.test(matrix(x,nr=2), alternative = "greater")$p.value)
  fisher_table$corr_p_value <- p.adjust(fisher_table$p_value, method = "BH")
  
  #Add a column with the ratio of selected sites to the ratio of total sites for each cluster
  fisher_table$SelectedToAllRatio <- fisher_table$SelectedNumInCluster/fisher_table$AllNumInCluster
  
  #Sort by corrected p-value
  fisher_table <- fisher_table[order(fisher_table$corr_p_value),]
  
  return(fisher_table)
}

#Takes output of cluster_enrichment function. Changes cluster name to be the first column rather than the row names.
#Selects just the cluster name and corrected p-value columns for clusters that meet a corrected p-value threshold.
#Output of this function can be used to join enrichment results of multiple experiments into a single table.
format_enrichment_table <- function(enrich_table, col_name, p_value) {
  #Cluster names are currently the row names. Make the cluster names the first column and delete the row names.
  #Name the new first column with the cluster names "Cluster"
  enrich_table <- cbind(rownames(enrich_table), data.frame(enrich_table, row.names=NULL))
  colnames(enrich_table)[1] <- "Cluster"
  
  #Select just the Cluster and corr_p_value columns.
  #Select the rows where the corr_p_value is less than the threshold p-value
  #Change name of corr_p_value column to col_name
  enrich_table <- enrich_table[c("Cluster", "corr_p_value")]
  enrich_table <- enrich_table[which (enrich_table$corr_p_value <  p_value), ]
  colnames(enrich_table)[2] <- col_name
  
  return(enrich_table)
}

#Takes a list of sites and a network file of CCCN edges. Plots site-site (i.e., CCCN) edges, site-gene edges (for those 
#sites where the parent gene is in the CFN), and gene-gene (i.e., CFN) edges (for the parent genes that are in the CFN)

graph.cluster <- function(sites, network) {
  cluster_edges <- filter.edges.0(sites, network)

  #Make gene-peptide edges
  peps <- unique(c(cluster_edges$source, cluster_edges$target))
  clust.cf <- gz.cf[which(gz.cf$id %in% peps),]
  gene_pep <- data.frame(source=clust.cf$Gene.Name, target=clust.cf$id, Weight=0.25, interaction="peptide")
  #Some gene nodes in the source column do not appear in the CFN; drop these rows
  gene_pep <- gene_pep[gene_pep$source %in% gz.cf$id,]

  #Get CFN edges between parent genes of cluster
  cfn_edges <- filter.edges.0(clust.cf$Gene.Name, gzalltgene.physical.cfn.merged)

  if (is.data.frame(cfn_edges)) {
    all.edges <- rbind(cluster_edges, gene_pep, cfn_edges)
    
  } else {
    all.edges <- rbind(cluster_edges, gene_pep)
  }

  #Get node attributes
  all.cf <- gz.cf[gz.cf$id %in% unique(c(all.edges$source, all.edges$target)),]
  
  #Plot the network
  cfn.cccn.suid <- createNetworkFromDataFrames(all.cf, all.edges, title=paste("CFN plus CCCN", (getNetworkCount()+1)), collection = "Interactions") 
  setNodeMapping(all.cf)
  setCorrEdgeAppearance(all.edges) 
  layoutNetwork("force-directed") 

}

#This function can be used in place of graph.cluster if a full CFN-CCCN has already been plotted in the Cytoscape session. 
#It skips the functions that set the node and edge properties (setNodeMapping and setCorrEdgeAppearance). 
#This makes drawing the networks much faster (10 seconds instead of 3 minutes). 
graph.cluster.fast <- function(sites, network) {
  cluster_edges <- filter.edges.0(sites, network)
  
  #Make gene-peptide edges
  peps <- unique(c(cluster_edges$source, cluster_edges$target))
  clust.cf <- gz.cf[which(gz.cf$id %in% peps),]
  gene_pep <- data.frame(source=clust.cf$Gene.Name, target=clust.cf$id, Weight=0.25, interaction="peptide")
  #Some gene nodes in the source column do not appear in the CFN; drop these rows
  gene_pep <- gene_pep[gene_pep$source %in% gz.cf$id,]
  
  #Get CFN edges between parent genes of cluster
  cfn_edges <- filter.edges.0(clust.cf$Gene.Name, gzalltgene.physical.cfn.merged)
  
  if (is.data.frame(cfn_edges)) {
    all.edges <- rbind(cluster_edges, gene_pep, cfn_edges)
    
  } else {
    all.edges <- rbind(cluster_edges, gene_pep)
  }
  
  #Get node attributes
  all.cf <- gz.cf[gz.cf$id %in% unique(c(all.edges$source, all.edges$target)),]
  
  #Plot the network
  cfn.cccn.suid <- createNetworkFromDataFrames(all.cf, all.edges, title=paste("CFN plus CCCN", (getNetworkCount()+1)), collection = "Interactions") 
  layoutNetwork("force-directed") 
  
}

#####################
#Takes a list of sites and a network file of CCCN edges. Creates Cytoscape edge and node tables that plot:
#site-site (i.e., CCCN) edges, site-gene edges (for those sites where the parent gene is in the CFN), 
#and gene-gene (i.e., CFN) edges (for the parent genes that are in the CFN). Writes the edge and node tables to files
#that can be manually loaded into Cytoscape. (Workaround for bug in RCy3 that was causing node and edge names to be permuted)
make.cluster.cytoscape.files <- function(sites, network, edge_filename, node_filename) {
  cluster_edges <- filter.edges.0(sites, network)
  
  #Make gene-peptide edges
  peps <- unique(c(cluster_edges$source, cluster_edges$target))
  clust.cf <- gz.cf[which(gz.cf$id %in% peps),]
  gene_pep <- data.frame(source=clust.cf$Gene.Name, target=clust.cf$id, Weight=0.25, interaction="peptide")
  #Some gene nodes in the source column do not appear in the CFN; drop these rows
  gene_pep <- gene_pep[gene_pep$source %in% gz.cf$id,]
  
  #Get CFN edges between parent genes of cluster
  cfn_edges <- filter.edges.0(clust.cf$Gene.Name, gzalltgene.physical.cfn.merged)
  
  if (is.data.frame(cfn_edges)) {
    all.edges <- rbind(cluster_edges, gene_pep, cfn_edges)
    
  } else {
    all.edges <- rbind(cluster_edges, gene_pep)
  }
  
  #Get node attributes
  all.cf <- gz.cf[gz.cf$id %in% unique(c(all.edges$source, all.edges$target)),]
  
  #Write the edge and node tables to files
  write.table(all.edges, file = edge_filename, sep = "\t", quote = FALSE, row.names=FALSE)
  write.table(all.cf, file = node_filename, sep = "\t", quote = FALSE, row.names=FALSE)
  
}

####################
#Takes a list of sites and a network file of CCCN edges. Gets the parent genes for all sites in the file.
#Creates Cytoscape edge and node tables that plots the CFN for these genes. Writes the edge and node tables to files
#that can be manually loaded into Cytoscape. (Workaround for bug in RCy3 that was causing node and edge names to be 
#permuted)

make.cfn.cytoscape.files <- function(sites, network, edge_filename, node_filename) {
  
  #Get CFN edges between parent genes of cluster
  sites.cf <- gz.cf[which(gz.cf$id %in% sites),]
  cfn_edges <- filter.edges.0(sites.cf$Gene.Name, gzalltgene.physical.cfn.merged)
  
  if (is.data.frame(cfn_edges)) {
    #Get node attributes
    genes.cf <- gz.cf[gz.cf$id %in% unique(c(cfn_edges$source, cfn_edges$target)),]
    
    #Write the edge and node tables to files
    write.table(cfn_edges, file = edge_filename, sep = "\t", quote = FALSE, row.names=FALSE)
    write.table(genes.cf, file = node_filename, sep = "\t", quote = FALSE, row.names=FALSE)
    
  } else {
    #If there are no cfn edges for the genes, just output the list of genes in a single column
    #Get the node attributes for these genes as well
    #This will be displayed as a set of unconnected nodes in Cytoscape
    gene_names_only <- data.frame(name = unique(sites.cf$Gene.Name))
      
    #Get node attributes
    genes.cf <- gz.cf[gz.cf$id %in% gene_names_only$name,]
    
    #Write the edge and node tables to files
    write.table(gene_names_only, file = edge_filename, sep = "\t", quote = FALSE, row.names=FALSE)
    write.table(genes.cf, file = node_filename, sep = "\t", quote = FALSE, row.names=FALSE)

  }
}

#####################
#For a given cluster, finds all the drug/cell-line combos where the cluster was significantly enriched 
#for drug changed sites and then returns the names of the replicate experiments for all of these drug/cell-line combos. 
#Takes a cluster id, a dataframe listing corrected p-values for each experiment that met a p-value threshold 
#(and NA for experiments that do not meet the threshold, e.g., enrichment_results_01.txt), and a list of the replicates
#for each experiment.
get_sig_experiment_reps <- function(cluster, enrichment_data, experiment_reps) {
  #Get the names of the significant experiments for the cluster
  selected_cluster <- enrichment_data[which(enrichment_data$Cluster == cluster),]
  selected_cluster[, c("Cluster", "num_exp")] <- list(NULL)
  selected_cluster_vector <- unlist(selected_cluster)
  sig_exp <- names(selected_cluster_vector[!is.na(selected_cluster_vector)])
  
  #Get the names of the replicates for all significant experiments
  selected_experiments <- unlist(experiment_reps[sig_exp])
  
  return(selected_experiments)
}


################
#Makes a table of PTM abundances for a given set of experiments and a given cluster appropriate for making a heatmape.
#Takes a cluster id, a dataframe listing the sites in each cluster (e.g., sites_by_cluster.txt), and a node data
#dataframe (e.g., gz.cf). Returns a dataframe with the sites in the cluster as rows and the abundance values 
#in all replicates of all significant experiments as columns. Abundance values of 0 are replaced with NA.
get_cluster_abundances <- function(cluster, selected_experiments, cluster_site_mapping, node_data) {
  
  #Select the sites (rows) that belong to the cluster from the node data table; select the PTM abundance columns
  #corresponding to the desired experiments from the node data table.
  cluster_sites <- cluster_site_mapping[ which(cluster_site_mapping$Cluster == cluster),]
  node_data_ptm <- node_data[which(node_data$Node.ID == "peptide"),] #Remove gene nodes from node data table
  node_data_ptm <- node_data[which(!(node_data$Gene.Name == "NA")),] #Remove nodes called NA from node data table
  
  node_data_ptm_cluster <- node_data_ptm[which(node_data_ptm$id %in% cluster_sites$Site),]
  cluster_data <- node_data_ptm_cluster[c("id", selected_experiments)]
  
  #Clean up the table:
  #To make the site names look nice, if there are multiple possible peptides for a site, keep only the first one
  #Make the site names the row names of the table
  #Replace the 0's with NA
  cluster_data$id <- sub(";.*", "", cluster_data$id) #If multiple possible PTMs, keep only the first one
  rownames(cluster_data) <- cluster_data$id
  cluster_data$id <- NULL
  cluster_data[cluster_data == 0] <- NA
  
  return(cluster_data)
}

##################
#Makes a heatmap where the rows and columns correspond to the rows and columns of the input data file (heatmap.data)
#Rows and columns are sorted by their means in descending orderl Color scheme is blue for low values and yellow for high
#Saves the heatmap as a 300dpi jpeg file
make.cluster.heatmap <- function(heatmap.data, file, image_width, image_height, order_cols = FALSE) {
  heatmap.data.m <- as.matrix(heatmap.data)
  #Order the rows and columns by the row means
  row_order <- order(-rowMeans(heatmap.data.m, na.rm = TRUE))
  if (order_cols == TRUE) {
    col_order <- order(-colMeans(heatmap.data.m, na.rm = TRUE))
    heatmap.data.mo <- heatmap.data.m[row_order, col_order, drop = FALSE]
  } else {
    heatmap.data.mo <- heatmap.data.m[row_order, ,drop = FALSE]
  }
  
  #Set a royal blue to yellow color palette
  rbyheatcolors <- colorRampPalette(colors=c('#0000FF',  '#FFFF00'), bias=0.25, space="rgb", interpolate = "spline")
  # royal blue to yellow
  palette(c(rbyheatcolors(500)))
  
  #Make heatmap and save to file
  jpeg(filename = file, units="in", width=image_width, height=image_height, res=300)
  if (ncol(heatmap.data.mo) == 1) {
    jpeg(filename = file, units="in", width=image_width, height=image_height, res=300)
    heatmap.2(cbind(heatmap.data.mo, heatmap.data.mo), dendrogram = "none", trace = "none", na.color="black", col=palette(),Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, symbreaks = TRUE)
    dev.off()
  } else {
    jpeg(filename = file, units="in", width=image_width, height=image_height, res=300)
    heatmap.2(heatmap.data.mo, dendrogram = "none", trace = "none", na.color="black", col=palette(),Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, symbreaks = TRUE)
    dev.off()
  }
  return(heatmap.data)
}

make.top.cluster.heatmap <- function(exp, cluster, cluster_site_mapping, node_data, filepath, filename, image_width, image_height) {
  reps <- unlist(experiments[exp])
  heatmap.data <- get_cluster_abundances(cluster, reps, cluster_site_mapping, node_data)
  file <- paste(filepath, filename, sep = "")
  make.cluster.heatmap(heatmap.data, file, image_width, image_height, order_cols = FALSE)
 
  
}


