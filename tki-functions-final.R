#################
#Calculate row medians for a dataframe; disregard NA values

rowMedians <- function(df) {
  
  meds <- apply(df, 1, median, na.rm = TRUE)
  return(meds)
}


#################
#If a vector has at least 2 (non-NA) values, return the median; otherwise, return NA.

rowMedians_w_reps <- function(vec) {
  if (sum(!is.na(vec)) > 1) {
    return(median(vec, na.rm = TRUE))
  } else {
    return(NA)
  }
}

##########
#Takes a dataframe as input. Returns the median of each row that has at least 2 (non-NA) values; returns NA for rows 
#with fewer than 2 values.

rowMedians_w_reps_df <- function(df) {
  meds <- apply(df, 1, rowMedians_w_reps)
  return(meds)
}

##################
#Selects sites whose value exceeds a threshold (up or down) in a particular condition
#Takes a dataframe of log2 values for different sites in different conditions, the name of the column(s) of values to use,
#a threshold value, and a filepath for saving the list of selected sites. Returns a dataframe containing the selected
#sites with a column for the site, a column for the value in the condition of interest and columns indicating 
#whether the site was up vs. down. The dataframe is also written to a file, and the number of selected sites is printed 
#to the console. Note: the input threshold is not a log2 value; log2 of the threshold is calculated within the function
#to be compatible with the log2 values in the dataframe. Note: when doing this analysis with medians of experiment groups
#site.df.exp will have only one column and count_up/count_down will be 0 (if median does not meet the threshold) or 1 (if median does
#meet the threshold)

get.sig.sites <- function(site.df, column, threshold = 2.25, filepath) {
  #Select the desired experiment from the sites table
  site.df.exp <- site.df[column]
  #Select sites that meet the threshold (up or down) for the desired experiment; rowSums is used to turn the logical into a numerical value
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

###############
#The cluster_enrichment function takes a dataframe showing which sites are in which clusters, a dataframe of 
#significantly changed sites, and a dataframe with a count of how many sites are in (and not in) each cluster. 
#It maps the sites to clusters, tallies the number of sig sites in each cluster (and not in each cluster), 
#calculates values for the 2x2 contingency matrix and performs the one-sided fisher test to determine whether there
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

######################
#The make.sig.site.heatmap.full function makes a heatmap of the significantly changed sites in a cluster.
#Significantly upregulated sites are yellow and significantly downregulated sites are blue. Inputs: cluster = the
#identifier of the cluster to be plotted; enrichment_data = a table of p-values for the clusters enriched in each experiment group
#sig_sites_list = a list where each element is a table of the significant sites for an experiment group; cluster_site_mapping =
#a two-column table showing the sites and the clusters they belong to; entire_cluster = if TRUE, all sites in the cluster
#will be plotted and if FALSE only the subset of sites provided in site_list will be plotted; all_exps = if TRUE all experiment
#groups will be plotted and if FALSE only those experiment groups with a significant p-value in enrichment_data will
#be plotted; filepath, filename = strings that specify the directory and filename to save the plot; image_height, image_width = 
#height and width of the heatmap

make.sig.site.heatmap.full <- function(cluster, enrichment_data, sig_sites_list, cluster_site_mapping, entire_cluster = TRUE, all_exps = FALSE, site_list = NULL, filepath, filename, image_width, image_height) {
  
  if (all_exps == FALSE) {
    #Only experiments groups for which the cluster is enriched will be plotted. Find all experiment groups where the cluster
    #is significantly enriched. Select that subset of experiment groups from the list of sig sites for all cell-line/drug combos
    sig.exps <- get_sig_experiments(cluster, enrichment_data)
    sig.sites.list <- sig_sites_list[sig.exps]
  } else {
    #All experiment groups plotted
    sig.sites.list <- sig_sites_list
  }
  
  #Classify the significant sites as up (1) or down (-1)
  sig.site.directions <- lapply(1:length(sig.sites.list), function(x) {get.directions(sig.sites.list[[x]])})
  names(sig.site.directions) <- names(sig.sites.list)
  
  #Select only the sig sites that are in the desired cluster
  sig.sites.in.cluster <- lapply(1:length(sig.site.directions), function(x) {join.with.cluster(sig.site.directions[[x]], cluster_site_mapping, cluster, names(sig.site.directions[x]))})
  names(sig.sites.in.cluster) <- names(sig.site.directions)
  
  #Convert the list (data for each experiment group is a separate list element) into a dataframe with one column
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
#The get_sig_experiments function takes the name of a cluster of interest (cluster) and a table showing enrichment p-values
#for each cluster for each experiment group (enrichment_data). Returns the names of the experiment groups for which the cluster
#has a significant p-value

get_sig_experiments <- function(cluster, enrichment_data) {
  #Get the names of the significant experiments for the cluster
  selected_cluster <- enrichment_data[which(enrichment_data$Cluster == cluster),]
  selected_cluster[, c("Cluster", "num_exp")] <- list(NULL)
  selected_cluster_vector <- unlist(selected_cluster)
  sig_exps <- names(selected_cluster_vector[!is.na(selected_cluster_vector)])
  
  return(sig_exps)
}

#################
#The get.directions function takes a table of sites that were significantly changed in an experiment group that has columns
#indicating the number of reps where the site was significantly up and the number of reps where the site was significantly down
#Returns the sites with a value of 1 if the number of sig up reps is greater than the number of sig down reps and
#a value of -1 otherwise. Also, if there are multiple names for a site separated by ";" only the first one is kept.
get.directions <- function(sig_sites_df) {
  sig_sites_df <- transform(sig_sites_df, direction = ifelse(count_up > count_down, 1 , -1))
  sig_sites_df <- sig_sites_df[c("id", "direction")]
  sig_sites_df$id <- sub(";.*", "", sig_sites_df$id)
  return(sig_sites_df)
}

###########################
#Takes a table of sites with a column indicating whether they are significantly up (1) or significantly down (-1),
#a table of sites indicating which cluster they are in, a cluster of interest, and an optional name for the column in the output
#indicating whether the site is up or down. Returns a table with a column showing the significant sites that are in the desired
#cluster and a column showing whether the site is significantly up (1) or down (-1).

join.with.cluster <- function(sig_site_directions, clusters, selected_cluster, col_name = "experiment") {
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
#Draws a heatmap with a blue (low) to yellow (high) diverging color palette and saves it to a file. 
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
  
  rbyheatcolors <- colorRampPalette(colors=c('#0000FF',  '#FFFF00'), bias=0.25, space="rgb", interpolate = "spline")
  # 100% blue to yellow
  palette(c(rbyheatcolors(500)))
  
  #Output file
  file <- paste(filepath, filename, sep = "")
  
  #Make heatmap and save to file
  
  if (ncol(heatmap.data.mo) == 1) {
    #This is an awkward workaround to handle cases where there is only one column of data for the heatmap because heatmap.2 won't
    #draw a single column heatmap. The single column is duplicated and the two identical columns are plotted side by side
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

##############
#The function make.cluster.heatmap.full makes a heatmap showing the abundance values for PTM sites in a cluster for desired
#experiments or groups of experiments. It can plot all sites in a cluster or a specified subset. It can plot individual experiments
#or aggregate by mean or median. It uses the function get_cluster_abundances to retrieve the abundances
#for all sites in a cluster for all selected experiments (or aggregates of experiments). Then it uses the function 
#make.cluster.heatmap to draw the heatmap.
#Inputs: cluster = cluster id; exp = a vector of names of experiments (or experiment groups) to be plotted; 
#rep_list = a list indicating the names of the reps for all experiment (names of reps should correspond to column names in 
#node_data); cluster_site_mapping = a dataframe listing the sites in each cluster (e.g., sites_by_cluster.txt)
#node_data = a dataframe that includes the abundances for all replicates of all experiments (e.g., gz.cf); aggregate_reps = if "mean"
#then the mean of the replicates for each experiment will be plotted, if "median" then the median of the replicates for each
#experiment will be plotted, and if "none" the individual replicates for each experiment will each be plotted in their own
#column; entire_cluster = if TRUE, all sites in cluster will be plotted and if FALSE, only sites specified in site_list will be
#plotted; filepath, filename = strings that specify the directory and filename to save the plot; image_height, image_width = 
#height and width of the heatmap

make.cluster.heatmap.full <- function(exp, rep_list, cluster, cluster_site_mapping, node_data, entire_cluster = TRUE, site_list = NULL, filepath, filename, image_width, image_height, aggregate_reps = "none", custom_order = FALSE, custom_order_vector = NULL, color = "yellow") {

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
#The function get_cluster_abundances returns a dataframe of PTM abundances for a given set of experiments and a given cluster 
#appropriate for making a heatmap. Inputs:
#cluster = cluster id; exp = a vector of names of experiments (or experiment groups) to be plotted; rep_list = a list indicating
#the names of the reps for all experiment (names of reps should correspond to column names in node_data); 
#cluster_site_mapping = a dataframe listing the sites in each cluster (e.g., sites_by_cluster.txt)
#node_data = a dataframe that includes the abundances for all replicates of all experiments (e.g., gz.cf); aggregate_reps = if "mean"
#then the mean of the replicates for each experiment will be plotted, if "median" then the median of the replicates for each
#experiment will be plotted, and if "none" the individual replicates for each experiment will each be plotted in their own
#column.

get_cluster_abundances <- function(cluster, exp, rep_list, cluster_site_mapping, node_data, aggregate_reps = "none") {
  
  #Select the sites (rows) that belong to the cluster from the node_data table; select the PTM abundance columns
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
  
  #Select the desired experiments from the list of all experiments and their reps
  selected_experiments <- rep_list[exp]
  
  if (aggregate_reps == "mean") {
    #Calculate the mean value of all reps for each experiment
    #Return a dataframe with the mean abundance over all reps for each site in the cluster for each experiment
    #selected_experiments is a list with experiment (e.g., h3122.criz) as name and character vector of names of replicates as
    #values. For each element in the list, calculate the rowMean of the columns from node_data_ptm_cluster whose names
    #correspond to the names of the replicates for that experiment
    cluster_data <- sapply(1:length(selected_experiments),function(x){rowMeans(node_data_ptm_cluster[selected_experiments[[x]]], na.rm = TRUE)})
    cluster_data <- as.data.frame(cluster_data)
    colnames(cluster_data) <- names(selected_experiments)
    rownames(cluster_data) <- rownames(node_data_ptm_cluster)
  } else if (aggregate_reps == "median") {
    #Calculate the median value of all reps for each experiment
    #Return a dataframe with the median abundance over all reps for each site in the cluster for each experiment
    #Similar to mean calculation above
    
    selected_experiments <- rep_list[exp]
    
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
#The function make.cluster.heatmap makes a heatmap where the rows and columns correspond to the rows and columns of the
#input data file and saves the heatmap as a 300dpi jpeg file.
#Inputs: heatmap.data = dataframe with values to be plotted; file = full path and name of file for saving heatmap; 
#image height/image width = height, width of plot; custom_order = if FALSE, rows will be plotted in descending order by row means
#and if FALSE, rows will be plotted in the order specified by custom_order_vector; color = if "red" plot will be colored using a 
#blue to red palette and if "yellow" plot will be colored using a blue to yellow palette

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
    # set a royal blue to yellow palette
    rbyheatcolors <- colorRampPalette(colors=c('#0000FF',  '#FFFF00'), bias=0.25, space="rgb", interpolate = "linear")
    palette(c(rbyheatcolors(500)))
  }
  
  #Make heatmap and save to file
  jpeg(filename = file, units="in", width=image_width, height=image_height, res=300)
  if (ncol(heatmap.data.mo) == 1) {
    #This is an awkward workaround to handle cases where there is only one column of data for the heatmap because heatmap.2 won't
    #draw a single column heatmap. The single column is duplicated and the two identical columns are plotted side by side
    heatmap.2(cbind(heatmap.data.mo, heatmap.data.mo), dendrogram = "none", trace = "none", na.color="black", col=rbyheatcolors,Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, cexRow = 1, symbreaks = TRUE)
  } else {
    heatmap.2(heatmap.data.mo, dendrogram = "none", trace = "none", na.color="black", col=rbyheatcolors,Rowv=NA, Colv=NA, hclustfun=NULL, density.info='none', margins=c(12, 10), srtCol = 45, offsetCol = 0, cexCol = 1, cexRow = 1, symbreaks = TRUE)
  }
  dev.off()
  return(heatmap.data)
}

####################
#The function plot_cluster_cfn plots the cfn connections of the parent genes for the sites in a cluster. Plotting is done 
#in Cytoscape with node shapes and borders and edge colors according to network style in PMID: 29789295. Inputs: cluster = name of cluster, 
#cluster_site_mapping = table mapping sites to clusters, ratio_table = table with metadata for each gene and site in the CCCN/CFN (like gz.cf or group_meds_all, which
#is like gz.cf except it provides median drug vs. vector ratios for several drug/cell-line experiment groups) and cfn = edge file
#of cfn interactions. 
plot_cluster_cfn <- function(cluster, cluster_site_mapping_table, ratio_table, cfn) {
  #select sites in cluster
  sites <- cluster_site_mapping_table[which (cluster_site_mapping_table$Cluster == cluster), "Site"]
  
  #Get parent gene names for sites in cluster
  genes <- ratio_table[which(ratio_table$id %in% sites), "Gene.Name"]
  
  #Get cfn edges for parent genes
  cfn_edges <- filter.edges.0(genes, cfn)
  
  #Get rows of ratio_table for each of the parent genes
  gene_info <- ratio_table[ratio_table$id %in% unique(c(cfn_edges$source, cfn_edges$target)),]
  
  #Plot the network in Cytoscape
  gene.suid <- createNetworkFromDataFrames(gene_info, cfn_edges, title=paste("Cluster", cluster, sep=" "), collection = "Interactions")
  setNodeMapping2(gene_info)
  edgeDprops.RCy32()
  layoutNetwork("force-directed")  
}

####################
#The function graph.cfn.ccn2 is like graph.cfn.cccn except it takes a file of node properties as input instead of always using
#gz.cf or ld.cf. Also simplifies id's with multiple genes to just show one gene

graph.cfn.cccn2 <- function(edgefile, nodefile = gz.cf.pruned, ld=FALSE, gz=TRUE, only.cfn=FALSE) {
  genenames <- extract.gene.names(edgefile)
  if (gz==TRUE) {
    cccn <- gzallt.cccnplus
    cccn$source <- sub(";.*", "", cccn$source)
    cccn$target <- sub(";.*", "", cccn$target)
    cccn.cf <- nodefile[nodefile$Gene.Name %in% genenames,]
  }
  if (ld==TRUE) {
    cccn <- ld.cccnplus
    cccn$source <- sub(";.*", "", cccn$source)
    cccn$target <- sub(";.*", "", cccn$target)
    cccn.cf <- nodefile[ld.cf$Gene.Name %in% genenames,]
  }
  if (only.cfn==TRUE) {
    cfn.cf <- cccn.cf[which(cccn.cf$Node.ID=="gene"),]
    gene.suid <- createNetworkFromDataFrames(cfn.cf, edgefile, title=paste("CFN", (getNetworkCount()+1)), collection = "Interactions")
    setNodeMapping2(cfn.cf)
    setCorrEdgeAppearance(edgefile)     
  }
  if (only.cfn==FALSE) {
    netpeps <- cccn.cf[which(cccn.cf$Node.ID=="peptide"), 'id']
    print(netpeps)
    # make gene-ptm edges
    net.gpe <- data.frame(source=cccn.cf$Gene.Name, target=cccn.cf$id, Weight=0.25, interaction="peptide")
    print(net.gpe)
    # remove gene-gene interactions (these aren't really autophosphorylations, just entries for the genes in gz.cf)
    net.gpe <- remove.autophos.RCy3(net.gpe)
    #make ptm-ptm edges
    ptm.cccn <-	filter.edges.0.RCy3(netpeps, cccn)
    print(ptn.cccn)
    cfn.cccn.edges <- rbind(net.gpe, ptm.cccn, edgefile)
    print(cfn.cccn.edges)
    if (gz==TRUE) {all.cf <- nodefile[nodefile$id  %in% unique(c(cfn.cccn.edges$source, cfn.cccn.edges$target)),]}
    if (ld==TRUE) {all.cf <- ld.cf[ld.cf$id  %in% unique(c(cfn.cccn.edges$source, cfn.cccn.edges$target)),]}
    print("nodes")
    print(all.cf$id)
    cfn.cccn.suid <- createNetworkFromDataFrames(all.cf, cfn.cccn.edges, title=paste("CFN plus CCCN", (getNetworkCount()+1)), collection = "Interactions") 
    print("start setNodeMapping2")
    setNodeMapping2(cccn.cf)
    print("start CorrEdgeAppearance")
    setCorrEdgeAppearance(cfn.cccn.edges) 
  }
  layoutNetwork("force-directed") 
  if (only.cfn==FALSE) return(cfn.cccn.edges)
}


####################
#The function plot_shortest_paths_cfn plots the shortest paths in the cfn between the parent genes of a set of sites and a 
#target gene. The function extracts the parent genes for the set of sites, calculates the shortest paths, extracts the
#metadata and ratio data for all genes in the shortest paths network and plots the network in Cytoscape. The node shape and border
#and edge styles are applied. (Coloring nodes according to ratio values happens separately.) The inputs are: drug_affected_sites =
#a dataframe with a column called "id" that lists the sites whose parent genes should be used for the shortest paths calculation;
#target = target gene(s) for the shortest path calculation; ratio_table = table with metadata for each gene and site in the CCCN/CFN (like gz.cf or group_meds_all, which
#is like gz.cf except it provides median drug vs. vector ratios for several drug/cell-line experiment groups), cfn = edge file
#of cfn interactions and title = string indicating the title for the network in Cytoscape. 
plot_shortest_paths_cfn <- function(drug_affected_sites, target, ratio_table, cfn, title) {
  #Get unique gene names for sig sites
  genes <- unique(ratio_table[which(ratio_table$id %in% drug_affected_sites$id), "Gene.Name"])
  
  #get shortest path edges
  shortest_path_edges <- composite.shortest.paths(genes, target, network = cfn)
  
  #Get rows of ratio_table for each of the parent genes
  gene_info <- ratio_table[ratio_table$id %in% unique(c(shortest_path_edges$source, shortest_path_edges$target)),]
  
  #Plot network; set style
  gene.suid <- createNetworkFromDataFrames(gene_info, shortest_path_edges, title=title, collection = "Shortest Paths")
  setNodeMapping2(gene_info)
  edgeDprops.RCy32()
  layoutNetwork("force-directed")  
}

###############################################
#Sets node shapes and border style; same as setNodeMapping except withe extra print statements

setNodeMapping2 <- function(cf) {
  print("start setBackgroundColorDefault")
  setBackgroundColorDefault(style.name = "default", new.color = "#949494") # grey 58
  print("start setNodeShapeDefault")
  setNodeShapeDefault("ELLIPSE")
  print("start setNodeColorDefault")
  setNodeColorDefault("#F0FFFF") # azure1
  print("start setNodeSizeDefault")
  setNodeSizeDefault(100) # for grey non-data nodes
  print("start setNodeFontSizeDefault")
  setNodeFontSizeDefault( 22)
  print("start setNodeLabelColorDefault")
  setNodeLabelColorDefault("#000000")  # black
  print("start setNodeBorderWidthDefault")
  setNodeBorderWidthDefault( 1.8)
  print("start setNodeBorderColorDefault")
  setNodeBorderColorDefault("#888888")  # gray 
  molclasses <- c("unknown", "receptor tyrosine kinase",  "SH2 protein", "SH2-SH3 protein", "SH3 protein", "tyrosine kinase",  "SRC-family kinase",   "kinase", "phosphatase", "transcription factor", "RNA binding protein")
  nodeshapes <- c("ELLIPSE","ROUND_RECTANGLE", "VEE", "VEE", "TRIANGLE", "HEXAGON", "DIAMOND", "OCTAGON", "OCTAGON", "PARALLELOGRAM", "RECTANGLE")
  print("start setNodeSelectionColorDefault")
  setNodeSelectionColorDefault("#CC00FF") 
  print("start setNodeShapeMapping")
  setNodeShapeMapping ("nodeType", molclasses, nodeshapes, default.shape="ELLIPSE")
  print("start setNodeBorderWidthMapping")
  setNodeBorderWidthMapping("nodeType", c("deacetylase","acetyltransferase","demethylase","methyltransferase","membrane protein", "receptor tyrosine kinase", "G protein-coupled receptor", "SRC-family kinase", "tyrosine kinase", "kinase", "phosphatase"), widths=c(4,12,4,12,8,16,16,12,12,12,14), 'd',default.width=4)
  if (length(cf[grep("SH2", cf$Domains), 1])>0 & !all(grep("SH2", cf$Domains) %in% which(cf$nodeType %in% molclasses))) {
    print("start setNodeShapeBypass SH2")
    setNodeShapeBypass(cf[grep("SH2", cf$Domains) %w/o% which(cf$nodeType %in% molclasses), 1], nodeshapes[3])} 
  if (length(cf[grep("RNA", cf$nodeType), 1])>0) {
    print("start setNodeShapeBypass RNA")
    setNodeShapeBypass(cf[grep("RNA", cf$nodeType), 1], nodeshapes[11])}
  if (length(cf[grep("transcription", cf$nodeType), 1])>0) {
    print("start setNodeShapeBypass transcription")
    setNodeShapeBypass(cf[grep("transcription", cf$nodeType), 1], nodeshapes[10])}
  if (length(cf[grep("acetyl", cf$nodeType), 1])>0) {
    print("start setNodeBorderColorBypass acetyl")
    setNodeBorderColorBypass(cf[grep("acetyl", cf$nodeType), 1], "#FF8C00")} # darkorange
  if (length(cf[grep("methyl", cf$nodeType), 1])>0) {
    print("start setNodeBorderColorBypass methyl")
    setNodeBorderColorBypass(cf[grep("methyl", cf$nodeType), 1], "#005CE6")} # blue
  if (length(cf[grep("membrane", cf$nodeType), 1])>0) {
    print("start set bypass membrane")
    print(cf[grep("membrane", cf$nodeType), 1])
    setNodeBorderColorBypass(cf[grep("membrane", cf$nodeType), 1], "#6600CC") # purple
    setNodeShapeBypass(cf[grep("membrane", cf$nodeType), 1], nodeshapes[2])} 
  if (length(cf[grep("kinase", cf$nodeType), 1])>0) {
    print("start setNodeBorderColorBypass kinase")
    setNodeBorderColorBypass(cf[grep("kinase", cf$nodeType), 1], "#EE0000")} # red2
  if (length(cf[grep("phosphatase", cf$nodeType), 1])>0) {
    print("start setNodeBorderColorBypass phosphatase")
    setNodeBorderColorBypass(cf[grep("phosphatase", cf$nodeType), 1], "#FFEC8B")} # lightgoldenrod1
  if (length(cf[grep("receptor", cf$nodeType), 1])>0) {
    print("start set bypass receptor")
    setNodeBorderColorBypass(cf[grep("receptor", cf$nodeType), 1], "#BF3EFF") # darkorchid1
    setNodeShapeBypass(cf[grep("receptor", cf$nodeType), 1], nodeshapes[2])} 
  if (length(cf[grep("TM", cf$nodeType), 1])>0) {
    print("start set bypass TM")
    setNodeBorderColorBypass(cf[grep("TM", cf$Domains), 1], "#6600CC") # purple
    setNodeShapeBypass(cf[grep("TM", cf$Domains), 1], nodeshapes[2])} 
}



############################
#%w/o% operator: return elements of x not in y

`%w/o%` <- function(x, y) {
  # return elements of x not in y
  x[ !(x %in% y) ]
}

################################
#Sets node fill color according to ratio values. Same as all.ratio.styles, but calls setNodeColorToRatios2
#which has appropriate size and color control points for log2 ratios
all.ratio.styles2 <- function(ratiocols=NULL) {
  nodevalues <- getTableColumns('node')
  
  if(length(ratiocols)==0) {
    ratiocolnames <- names(nodevalues)[grep("atio", names(nodevalues))] %w/o% "No.Modifications"} else {ratiocolnames <- names(nodevalues[ratiocols])}

  for (i in 1:length(ratiocolnames)){ 
    plotcol <- ratiocolnames[i]
    print(paste("Plotcol: ", plotcol, sep = ""))
    style.name = paste(ratiocolnames[i], "Style")
    print(style.name)
    setVisualStyle("default")
    setNodeColorToRatios2(plotcol)    
    copyVisualStyle('default', style.name)
    setVisualStyle(style.name)
  }
} 

################################
#Sets node size and color to match ratio data using columns from Cytoscape node table. Called by all.ratio.styles
#Same as setNodeColorToRatios except size and color control points are appropriate for log2 ratios
setNodeColorToRatios2 <- function(plotcol){
  cf <- getTableColumns('node')
 
  if(!(plotcol %in% getTableColumnNames('node'))){
    
    cat("\n","\n","\t", "Which attribute will set node size and color?")
    plotcol <- as.character(readLines(con = stdin(), n = 1))
  }
  
  node.sizes = c(135, 130, 108, 75, 35, 75, 108, 130, 135)
  ratio.colors = c ('#0099FF', '#007FFF','#00BFFF', '#00CCFF', '#00FFFF', '#00EE00', '#FFFF7E', '#FFFF00', '#FFE600', '#FFD700', '#FFCC00')
  #size control points for original ratios were: -100.0, -15.0, -5.0, 0.0, 5.0, 15.0, 100.0, but we are using log2 ratios
  size.control.points = c(-log(100.0, 2), -log(15.0, 2), -log(5.0,2), 0.0, log(5.0, 2), log(15.0, 2), log(100.0, 2))
  
  #color control points for original ratios were: -100.0, -10.0, -5.0, -2.25, 0.0, 2.25, 5.0, 10.0, 100.0, but we are using log2 ratios
  color.control.points = c(-log(100.0, 2), -log(10.0, 2), -log(5.0, 2), -log(2.25, 2), 0.0, log(2.25, 2), log(5.0, 2), log(10.0, 2), log(100.0, 2))
  
  limits <- range(cf[, plotcol], na.rm = TRUE)
  if(limits[1] < min(size.control.points)) {
    size.control.points[1] = limits[1] - 0.1
  }
  
  if(limits[1] < min(color.control.points)) {
    color.control.points[1] = limits[1] - 0.1
  }
  
  if(limits[2] > max(size.control.points)) {
    size.control.points[length(size.control.points)] = limits[2] + 0.1
  }
  
  if(limits[2] > max(color.control.points)) {
    color.control.points[length(color.control.points)] = limits[2] + 0.1
  }

  setNodeColorMapping(names(cf[plotcol]), color.control.points, ratio.colors, 'c')
  lockNodeDimensions('TRUE')
  setNodeSizeMapping (names(cf[plotcol]), size.control.points, node.sizes, 'c')
  setNodeSelectionColorDefault ( "#CC00FF") 
}
