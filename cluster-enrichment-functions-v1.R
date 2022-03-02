###############
#Functions for Cluster Enrichment Analysis
##############

#############################
#The function get.sig.sites takes a list of columns from gz.cf (assumes that the column list is "id", "Gene.Name", 
#and then the ratio columns for the desired drug/cell line combo), a threshold fold-change, and a minimum experiment
#number and returns a dataframe with the id, gene name, and values for sites that are changed by at least the threshold
#in the same direction in at least the minimum specified number of experiments. If the total number of ratio columns
#is less than the specified minimum, then it will return sites that meet the threshold for all of the given ratio 
#columns. For example, if num_exp = 3, but there are only two ratio colums given, the function will return sites 
#which are changed by at least the threshold in both of the ratio columns.

get.sig.sites <- function(site_df = gz.cf, columns, threshold = 2.25, num_exp = 2) {
  #Remove gene lines from gz.cf
  gz.cf.sites <- gz.cf[which (gz.cf$Node.ID == "peptide"), ]
  #Make a smaller version of gz.cf with only the site, gene, and desired ration columns
  gz.cf.selected <- gz.cf.sites[columns]
  #Count sites that meet the threshold (up or down)
  gz.cf.selected$count_up <- rowSums(gz.cf.selected[3:length(columns)] > log2(threshold))
  gz.cf.selected$count_down <- rowSums(gz.cf.selected[3:length(columns)] < -log2(threshold))
  num_ratio_cols <- length(columns) - 2
  if (num_ratio_cols < num_exp) {
    #Require that threshold be met for all ratio columns
    selected.sig.sites <- gz.cf.selected[which((gz.cf.selected$count_down == num_ratio_cols) | (gz.cf.selected$count_up == num_ratio_cols)),]
  } else {
    #Require that the threshold be met for at least the number of ratio columns specified in num_exp
    selected.sig.sites <- gz.cf.selected[which((gz.cf.selected$count_down > (num_exp - 1)) | (gz.cf.selected$count_up > (num_exp - 1))),]
  }
  
  return(selected.sig.sites) }


###############
#The cluster_enrichment function takes a dataframe showing which sites are in which clusters, a dataframe of 
#significantly changed sites, and a dataframe with a count of how many sites are in (and not in) each cluster. 
#It maps the sites to clusters, tallies the number of sig sites in each cluster (and not in each cluster), 
#cacluates values for the 2x2 contingency matrix and performs the one-sided fisher test to determine whether there
#are more significantly changed sites in a cluster than expected by chance. Returns a dataframe with the values used 
#for the contingency table, the raw p-value, the Benjamini corrected p-value and the fraction of sites in the cluster
#that are significantly changed.

cluster_enrichment <- function(clusters, sig_sites, all_sites_tally) {
  #Change name of id column of sig_sites table to "Site" so it can be merged with the clusters table
  names(sig_sites)[names(sig_sites) == 'id'] <- 'Site'
  
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
  
  #Merge the data for all clusters and the data for the selected clusters. Calculate the other two values needed for
  #the 2x2 contingency table: the non-significant sites that are in the cluster (NotSelectedNumInCluster) and
  #the non-significant sites that are not in the cluster (NotSelectedNumOutCluster). Remove the two all cluster
  #columns (AllNumInCluster and  AllNumOutCluster) so only the 4 columns needed for the contingency table remain.
  #Set the row names to be the cluster ID's. Remove rows where the number of selected sites in a cluster
  #is less than 3 because this is unlikely to give a significant result and will just drive up the number of tests
  fisher_table <- merge(selected_tally, all_sites_tally, by = "Cluster")
  fisher_table$NotSelectedNumInCluster <- fisher_table$AllNumInCluster - fisher_table$SelectedNumInCluster
  fisher_table$NotSelectedNumOutCluster <- fisher_table$AllNumOutCluster - fisher_table$SelectedNumOutCluster
  fisher_table$AllNumInCluster <- NULL
  fisher_table$AllNumOutCluster <- NULL
  rownames(fisher_table) <- fisher_table$Cluster
  fisher_table$Cluster <- NULL
  fisher_table <- fisher_table[ which(fisher_table$SelectedNumInCluster >= 3),]
  
  
  #Perform the Fisher Test and the Benjamini-Hochberg multiple testing correction
  fisher_table$p_value <- apply(fisher_table,1, function(x) fisher.test(matrix(x,nr=2), alternative = "greater")$p.value)
  fisher_table$corr_p_value <- p.adjust(fisher_table$p_value, method = "BH")
  
  #Add a column with the ratio of selected sites to the ratio of total sites for each cluster
  fisher_table$SelectedToAllRatio <- fisher_table$SelectedNumInCluster/(fisher_table$SelectedNumInCluster + fisher_table$NotSelectedNumInCluster)
  #Sort by corrected p-value
  fisher_table <- fisher_table[order(fisher_table$corr_p_value),]
  
  #fisher_table$Cluster <- rownames(fisher_table)
  #rownames(fisher_table) <- NULL
  return(fisher_table)
}

######################
#The format_enrichment_table function takes the dataframe output by the cluster_enrichment function and returns
#a two-column dataframe with the cluster number in the first column and the corrected p-values in the second column.
#Only rows (clusters) that meet the desired corrected p-value threshold are retained.

format_enrichment_table <- function(enrich_table, col_name, p_value) {
  #Cluster names are currently the row names. Make the cluster names the first column and deletes the row names.
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

#######################
# The function filter.edges.0 takes a list of node names and an edge file and selects edges from the edge file
#that contain at least one of the nodes on the input list

filter.edges.0 <- function(nodenames, edge.file) {
  nodenames <-as.character(nodenames)
  a = as.character(edge.file[,1])
  b = as.character(edge.file[,2])
  edgefile.nodes <- unique(c(a,b))
  sel.edges <- edge.file[edge.file[,1] %in% nodenames & edge.file[,2] %in% nodenames,]
  if(dim(sel.edges)[1] == 0) {return(NA)} else return(sel.edges) 
}

####################
#The function make.cfn.cytoscape.files takes a list of sites. Gets the parent genes
#for all sites in the file. Creates Cytoscape edge and node tables that plots the CFN for these genes. 
#Writes the edge and node tables to files that can be manually loaded into Cytoscape. (Workaround for bug in RCy3 that
#was causing node and edge names to be permuted)

make.cfn.cytoscape.files <- function(sites, edge_filename, node_filename) {
  
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

##############
#The function make.top.cluster.heatmap takes name(s) of one or more drug/cell-line combos. It gets the names of the
#replicates for the drug/cell-line combos. Then it uses the function get_cluster_abundances to retrieve the abundances
#for all sites in a cluster for all selected replicates. Then it uses the function make.cluster.heatmap to draw a 
#heatmap of the abundances for all selected replicates.
make.top.cluster.heatmap <- function(exp, rep_list, cluster, cluster_site_mapping, node_data, entire_cluster = TRUE, site_list = NULL, filepath, filename, image_width, image_height, aggregate_reps = "none", custom_order = FALSE, custom_order_vector = NULL, color = "yellow") {
  reps <- unlist(rep_list[exp])
  heatmap.data <- get_cluster_abundances(cluster, exp, rep_list, cluster_site_mapping, node_data, aggregate_reps)
  
  #Optional: If entire_cluster = FALSE and a vector of sites is provided, select the subset of rows from the data table
  #where the site id (row name) is on the selected site list
  if (entire_cluster == FALSE) {
      heatmap.data <- heatmap.data[ which(rownames(heatmap.data) %in% site_list),]
  }
  
  file <- paste(filepath, filename, sep = "")
  make.cluster.heatmap(heatmap.data, file, image_width, image_height, custom_order, custom_order_vector, color = color)
  
  
}

################
#Makes a table of PTM abundances for a given set of experiments and a given cluster appropriate for making a heatmape.
#Takes a cluster id, a dataframe listing the sites in each cluster (e.g., sites_by_cluster.txt), and a node data
#dataframe (e.g., gz.cf). If mean_values is FASLE, returns a dataframe with the sites in the cluster as rows and 
#the abundance values in all replicates of all significant experiments as columns. If mean_values = TRUE, the mean
#of replicates is returned for each experiment. Abundance values of 0 are replaced with NA.
get_cluster_abundances <- function(cluster, exp, rep_list, cluster_site_mapping, node_data, aggregate_reps = "none") {
  
  #Select the sites (rows) that belong to the cluster from the node data table; select the PTM abundance columns
  #corresponding to the desired experiments from the node data table.
  cluster_sites <- cluster_site_mapping[ which(cluster_site_mapping$Cluster == cluster),]
  node_data_ptm <- node_data[which(node_data$Node.ID == "peptide"),] #Remove gene nodes from node data table
  node_data_ptm <- node_data[which(!(node_data$Gene.Name == "NA")),] #Remove nodes called NA from node data table
  
  node_data_ptm_cluster <- node_data_ptm[which(node_data_ptm$id %in% cluster_sites$Site),]
  
  #Clean up the table:
  #To make the site names look nice, if there are multiple possible peptides for a site, keep only the first one
  #Make the site names the row names of the table
  #Replace the 0's with NA
  node_data_ptm_cluster$id <- sub(";.*", "", node_data_ptm_cluster$id) #If multiple possible PTMs, keep only the first one
  rownames(node_data_ptm_cluster) <- node_data_ptm_cluster$id
  node_data_ptm_cluster$id <- NULL
  node_data_ptm_cluster[node_data_ptm_cluster == 0] <- NA
  
  if (aggregate_reps == "mean") {
    #Calculate the mean value of all reps for each experiment
    #Return a dataframe with the mean abundance over all reps for each site in the cluster for each experiment
    
    selected_experiments <- rep_list[exp]
    
    #selected_experiments is a list with experiment (e.g., h3122.criz) as name and character vector of names of replicates as value
    #for each element in the list, calculate the rowMean of the columns from node_data_ptm_cluster whose names
    #correspond to the names of the replicates for that experiment
    cluster_data <- sapply(1:length(selected_experiments),function(x){rowMeans(node_data_ptm_cluster[selected_experiments[[x]]], na.rm = TRUE)})
    cluster_data <- as.data.frame(cluster_data)
    colnames(cluster_data) <- names(selected_experiments)
    rownames(cluster_data) <- rownames(node_data_ptm_cluster)
  } else if (aggregate_reps == "median") {
    #Calculate the median value of all reps for each experiment
    #Return a dataframe with the median abundance over all reps for each site in the cluster for each experiment
    
    selected_experiments <- rep_list[exp]
    
    #selected_experiments is a list with experiment (e.g., h3122.criz) as name and character vector of names of replicates as value
    #for each element in the list, calculate the row medians of the columns from node_data_ptm_cluster whose names
    #correspond to the names of the replicates for that experiment
    cluster_data <- sapply(1:length(selected_experiments),function(x){rowMedians(node_data_ptm_cluster[selected_experiments[[x]]])})
    cluster_data <- as.data.frame(cluster_data)
    colnames(cluster_data) <- names(selected_experiments)
    rownames(cluster_data) <- rownames(node_data_ptm_cluster)
  } else {
    
    #Return a dataframe with abundance for each individual replicate for each site in the cluster for each experiment
   
    reps <- unlist(rep_list[exp])
    cluster_data <- node_data_ptm_cluster[reps]
  }

  
  return(cluster_data)
}

##################
#The function make.cluster.heatmap makes a heatmap where the rows and columns correspond to the rows and columns of the input data file (heatmap.data)
#Rows and columns are sorted by their means in descending order. Color scheme is blue for low values and red for high
#Saves the heatmap as a 300dpi jpeg file
make.cluster.heatmap <- function(heatmap.data, file, image_width, image_height, custom_order = FALSE, custom_order_vector = NULL, color = "yellow") {
  
  
  if (custom_order == FALSE) {
    #Order the rows by the row means
    row_order <- order(-rowMeans(heatmap.data, na.rm = TRUE))
    heatmap.data.o <- heatmap.data[row_order, ,drop = FALSE]
  } else {
    #Order the rows using a custom vector
    heatmap.data.o <- heatmap.data %>% arrange(factor(rownames(heatmap.data), levels = custom_order_vector))
  }
  
  heatmap.data.mo <- as.matrix(heatmap.data.o)
  
  if (color == "red") {
  #Set a blue to red diverging color palette
  rbyheatcolors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255)
  } else {
    # set a royal blue to palette
    rbyheatcolors <- colorRampPalette(colors=c('#0000FF',  '#FFFF00'), bias=0.25, space="rgb", interpolate = "linear")
    palette(c(rbyheatcolors(500)))
  }
  
  #Make heatmap and save to file
  jpeg(filename = file, units="in", width=image_width, height=image_height, res=300)
  if (ncol(heatmap.data.mo) == 1) {
    jpeg(filename = file, units="in", width=image_width, height=image_height, res=300)
    #heatmap.2(cbind(heatmap.data.mo, heatmap.data.mo), dendrogram = "none", trace = "none", na.color="black", col=viridis(15),Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, cexRow = 1, symbreaks = TRUE)
    heatmap.2(cbind(heatmap.data.mo, heatmap.data.mo), dendrogram = "none", trace = "none", na.color="black", col=rbyheatcolors,Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, cexRow = 1, symbreaks = TRUE)
    dev.off()
  } else {
    jpeg(filename = file, units="in", width=image_width, height=image_height, res=300)
    #heatmap.2(heatmap.data.mo, dendrogram = "none", trace = "none", na.color="black", col=viridis(15),Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, cexRow = 1, symbreaks = TRUE)
    heatmap.2(heatmap.data.mo, dendrogram = "none", trace = "none", na.color="black", col=rbyheatcolors,Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, cexRow = 1, symbreaks = TRUE)
    dev.off()
  }
  return(heatmap.data)
}

######################
#Makes a heatmap of the significantly changed sites in a cluster for all cell-line/drug combos in which the cluster
#is enriched for sig changed sites. Significantly upregulated sites are red and significantly downregulated sites
#are blue. Inputs are: the cluster of interest, a table of p-values for the clusters enriched in each cell-line/drug
#combo, a list where each element is a table of the significant sites for a cell-line/drug combo, the desired file pathand name for the output heatmap file, and the desired width and height of the heatmap.


#make.sig.site.heatmap.full("38.41.40", enrich.table, sig.site.meds.list, CCCN_sites, all_exps = TRUE
make.sig.site.heatmap.full <- function(cluster, enrichment_data, sig_sites_list, cluster_sites, entire_cluster = TRUE, all_exps = FALSE, site_list = NULL, filepath, filename, image_width, image_height) {
  
  if (all_exps == FALSE) {
    #Find all cell-line/drug combos where the cluster is significantly enriched. Select that subset of cell-line/drug
    #combos from the list of sig sites for all cell-line/drug combos
    sig.exps <- get_sig_experiments(cluster, enrichment_data)
    sig.sites.list <- sig_sites_list[sig.exps]
  } else {
    sig.sites.list <- sig_sites_list
  }
  
  #Classify the sites as significantly up (1) or significantly down (-1)
  sig.site.directions <- lapply(1:length(sig.sites.list), function(x) {get.directions(sig.sites.list[[x]])})
  names(sig.site.directions) <- names(sig.sites.list)
  
  #Select only the sig sites that are in the desired cluster
  sig.sites.in.cluster <- lapply(1:length(sig.site.directions), function(x) {join.with.cluster(sig.site.directions[[x]], cluster_sites, cluster, names(sig.site.directions[x]))})
  names(sig.sites.in.cluster) <- names(sig.site.directions)
  
  #Convert the list (data for each cell-line/drug combo is a separate list element) into a dataframe with one column
  #for each cell-line/drug combo 
  sig.sites.df <- Reduce(merge, sig.sites.in.cluster)
  
  #Optional: only keep a subset of the sites in a cluster. In this case entire_cluster = FALSE, and a site list should
  #be provided 
  if (entire_cluster == FALSE) {
    sig.sites.df <- pick.selected.sites(site_list, sig.sites.df)
  }
  
  #make the site id's into the row names
  rownames(sig.sites.df) <- sig.sites.df$id
  sig.sites.df$id <- NULL
  
  #Make the heatmap
  heatmap.final <- make.sig.site.heatmap(sig.sites.df, filepath, filename, image_width, image_height)
  return(heatmap.final)
  
}

###################
#Takes the name of a cluster of interest and a table showing enrichment p-values for each cluster for each
#cell-line/drug combo. Returns the names of the cell-line/drug combos for which the cluster has a significant p-value
get_sig_experiments <- function(cluster, enrichment_data) {
  #Get the names of the significant experiments for the cluster
  selected_cluster <- enrichment_data[which(enrichment_data$Cluster == cluster),]
  selected_cluster[, c("Cluster", "num_exp")] <- list(NULL)
  selected_cluster_vector <- unlist(selected_cluster)
  sig_exps <- names(selected_cluster_vector[!is.na(selected_cluster_vector)])
  
  return(sig_exps)
}

#################
#Takes a table of sites that were significantly changed in a cell-line/drug combos that has columns indicating
#the number of reps where the site was significantly up and the number of reps where the site was significantly down
#Returns the sites with a value of 1 if the number of sig up reps is greater than the number of sig down reps and
#a value of -1 otherwise.
get.directions <- function(sig_sites_df) {
  sig_sites_df <- transform(sig_sites_df, direction = ifelse(count_up > count_down, 1 , -1))
  sig_sites_df <- sig_sites_df[c("id", "direction")]
  sig_sites_df$id <- sub(";.*", "", sig_sites_df$id)
  return(sig_sites_df)
}


###########################
#Takes a table of sites with a column indicating whether they are significantly up (1) or significantly down (-1),
#a table of sites indicating which cluster they are in, and a cluster of interest. Returns the subset of rows where
#the sites are in the desired cluster.

join.with.cluster <- function(sig_site_directions, clusters, selected_cluster, col_name = "test3") {
  clusters <- rename(clusters, id = Site)
  clusters$id <- sub(";.*", "", clusters$id)
  cluster_sites <- full_join(sig_site_directions, clusters, by = "id")
  cluster_sites <- cluster_sites[ which(cluster_sites$Cluster %in% selected_cluster),c("id", "direction")]
  colnames(cluster_sites)[2] <- col_name
  return(cluster_sites)
  
}

#######################
#Takes a dataframe of with a column of site ids (called id) and one or more columns of values for the sites
#in different experiments and a vector of site ids of interest. Returns the rows of the df where the site
#is in the list of interest
pick.selected.sites <- function(site_list, all_sites) {
  selected_sites <- all_sites[ which(all_sites$id %in% site_list),]
  return(selected_sites)
}

###################
#Draws a heatmap with a blue to red diverging color palette and saves it to a file. 
#Orders the heatmap rows according to the sum of the values in the row with rows that are all NA at the bottom.
make.sig.site.heatmap <- function(heatmap.data, filepath, filename, image_width, image_height) {
  heatmap.data.m <- as.matrix(heatmap.data)
  #Order the rows by the row means multiplied by the number of elements in the row that are not NA; this is equivalent
  #to doing ordering by row sums except for the way NAs are handled. If you just do rowSums, rows that are all NA will
  #have a value of 0; if you do row means multiplied by number of elements, rows that are all NA have a value of NA
  #(and will get sorted to the bottom of the plot)
  sums_na <- rowMeans(heatmap.data.m, na.rm = TRUE) * rowSums(!is.na(heatmap.data.m))
  row_order <- order(-sums_na)
  heatmap.data.mo <- heatmap.data.m[row_order, ,drop = FALSE]
  
  #Set a blue to red diverging color palette
  #rbyheatcolors <- colorRampPalette( rev(brewer.pal(5, "RdBu")) )(255)
  
  rbyheatcolors <- colorRampPalette(colors=c('#0000FF',  '#FFFF00'), bias=0.25, space="rgb", interpolate = "spline")
  # 100% blue to yellow
  palette(c(rbyheatcolors(500)))
  
  #Output file
  file <- paste(filepath, filename, sep = "")
  
  #Make heatmap and save to file
  
  if (ncol(heatmap.data.mo) == 1) {
    jpeg(filename = file, units="in", width=image_width, height=image_height, res=300)
    heatmap.2(cbind(heatmap.data.mo, heatmap.data.mo), dendrogram = "none", trace = "none", na.color="black", col=rbyheatcolors,Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, symbreaks = TRUE)
    dev.off()
  } else {
    jpeg(filename = file, units="in", width=image_width, height=image_height, res=300)
    heatmap.2(heatmap.data.mo, dendrogram = "none", trace = "none", na.color="black", col=rbyheatcolors,Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, cexRow = 1, symbreaks = TRUE)
    dev.off()
  }
  
  heatmap.data.o.df <- as.data.frame(heatmap.data.mo)
  return(heatmap.data.o.df)
}

#################
rowMedians <- function(df) {
   
    meds <- apply(df, 1, median, na.rm = TRUE)
    return(meds)
}

#################
rowMedians_w_reps <- function(vec) {
    if (sum(!is.na(vec)) > 1) {
      return(median(vec, na.rm = TRUE))
    } else {
        return(NA)
    }
}

#########
rowMedians_w_reps_df <- function(df) {
  meds <- apply(df, 1, rowMedians_w_reps)
  return(meds)
}

##################
get.sig.sites.meds <- function(site.df = group_meds, column, threshold = 2.25, filepath = sig.site.meds.filepath) {
  #Select the desired experiment from the sites table
  site.df.exp <- site.df[column]
  #Select sites that meet the threshold (up or down) for the desired experiment
  site.df.exp$count_up <- rowSums(site.df.exp[1] > log2(threshold))
  site.df.exp$count_down <- rowSums(site.df.exp[1] < -log2(threshold))
  site.df.sig <- site.df.exp[which((site.df.exp$count_down == 1) | (site.df.exp$count_up == 1)),]
  
  #Print the number of sig sites
  print(paste(column, ": ", nrow(site.df.sig), sep = ""))
  
  #Make a column called 'id' with the sites (currently the rownames); this is to make the input compatible with the
  #cluster enrichment function
  site.df.sig$id <- rownames(site.df.sig)
  rownames(site.df.sig) <- NULL
  out.file <- paste(filepath, column, ".sig.sites.meds.txt", sep = "")
  write.table(site.df.sig, file = out.file, sep = "\t", quote = FALSE, col.names=NA)
  return(site.df.sig)
}