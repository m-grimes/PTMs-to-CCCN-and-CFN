# Protein-Protein Interaction Networks
# July 5, 2019
# Mark Grimes

# =============================================================================
# The next step is to use the gene CCCN to filter PPI edges for a PPI CFN
#______________________________________________________________________________
# =======================     
# Retrive PPI edges from several sources.
#_______________________________________________________________
#
# load work
load(file=paste("_LINCS/_KarenGuolin/", "GZ_PPI_Networks2.RData", sep=""))
# load("/Users/_mark_/Dropbox/_Work/R_/_LINCS/Bin_Curated_Data/ld_lim_PMD.RData")
load(file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "TenCell.RData", sep=""))
#_______________________________________________________________
# Update Pathway Commons database
# Pathway commons retrieved from http://www.pathwaycommons.org/archives/PC2/v9/
# pcname <- "PathwayCommons9.All.hgnc.txt"
# New version 11 https://www.pathwaycommons.org/archives/PC2/v11/
pcname <- "PathwayCommons11.All.hgnc.txt"
pcpath <- "/Volumes/Terra_Byte/R_Archive_2/PC_BioPlex/"
pcpath <- "/Users/mark_grimes/Software\ Images/" # Laptop
# load v9 version
# load(file=paste(pcpath,"pcommons_2018.RData", sep=""))
load(file=paste(pcpath,"pcommons_2019.RData", sep=""))
pcfile <- paste(pcpath, pcname, sep="")
pcommons <- read.table(pcfile, header=TRUE, sep = "\t", comment.char = "#", na.strings='', stringsAsFactors=FALSE, fill=TRUE)
# Dim  1014008       7
unique(pcommons$INTERACTION_DATA_SOURCE[grep("PhosphoSite", pcommons$INTERACTION_DATA_SOURCE)])
# Note: Now Contains PSP info
unique(pcommons$INTERACTION_DATA_SOURCE[grep("lex", pcommons$INTERACTION_DATA_SOURCE)])
# But does not contain BioPlex info.
# Now contains interactions with DNA, RNA, and transcription complexes
names(pcommons)
unique(pcommons$INTERACTION_TYPE)
head(pcommons[grep("DnaReference", pcommons$INTERACTION_TYPE),1:6])
head(pcommons[grep("RnaReference", pcommons$INTERACTION_TYPE),1:6])
head(pcommons[grep("ProteinReference", pcommons$INTERACTION_TYPE),1:6])
# Note the format is different 
head(pcommons[grep(";", pcommons$PARTICIPANT_B),1:6])
# But the DNA and RNA interactions can be excluded using 
head(pcommons[-grep(";", pcommons$PARTICIPANT_B),1:6])
# There are also metabolites and chemical compounds
head(pcommons[grep("reacts-with" , pcommons$INTERACTION_TYPE),1:6])
head(pcommons[grep("used-to-produce" , pcommons$INTERACTION_TYPE),1:6])
head(pcommons[grep("consumption-controlled-by" , pcommons$INTERACTION_TYPE),1:6])
head(pcommons[grep("chemical-affects" , pcommons$INTERACTION_TYPE),1:6])
head(pcommons[grep("controls-production-of", pcommons$INTERACTION_TYPE),1:6])
# These have names in the format CHEBI:NNNN for compounds.
# There is no longer "neighbor-of"
head(pcommons[grep("neighbor" , pcommons$INTERACTION_TYPE),1:6])
# Convert to Formatting of PPI edge files used in my functions is
# "Gene.1"   "Gene.2"   "Weight"   "edgeType"
# But note that Cytoscape now likes the edge names to be
# "source"     "target"  "Weight"  "interaction"
# So this will have to be converted for graphing in Cytoscape
pcnet <- pcommons[-grep("CHEBI:", pcommons$PARTICIPANT_A),1:4]
pcnet <- pcnet[-grep("CHEBI:", pcnet$PARTICIPANT_B),]
# Add an arbitrary Weight and rearrange to be consistent with other PPI networks
pcnet$Weight <- 0.2
pcnet <- pcnet[, c(1,3,5,2,4)]
names(pcnet) <- c("Gene.1", "Gene.2", "Weight", "edgeType", "Source")
# now dim 633126      5 (v9);  594752      5 (v11) 
# NOTE: still includes gene and mRNA interactions, which may be useful later
# save(pcnet, file=paste(comp_path,"Dropbox/_Work/R_/pcnet.RData", sep=""))
# ________________________________#________			
#  New Network
#  		BioPlex: Huttlin et al., 2015, Cell 162, 425–440
#		parameters: incorrect ID (p(Wrong)=pW), background (p(NoInt)=pNI), or a true interactor (p(Int)=pInt) 

bpfilename <- "/Volumes/Terra_Byte/R_Archive_2/PC_BioPlex/BioPlex_interactionList_v4.tsv"
bioplex <- read.table(bpfilename, header=TRUE, sep = "\t", comment.char = "", na.strings='', stringsAsFactors=FALSE, fill=TRUE)
# 56553 edges; pInt values range from 0.75 to 1.00
bioplex <- bioplex[,c("SymbolA", "SymbolB", "pInt")]
names(bioplex) <- c("Gene.1", "Gene.2", "Weight")
bioplex$edgeType <- "BioPlex"
save(pcommons, bioplex, pcnet, enzsub, file=paste(pcpath,"pcommons_2019.RData", sep=""))
#_______________________________________________________________
# Karen's enzyme substrate data
enzsubname <- "iptmnet4.1_enzyme_substrate_relations_scored_simpler.txt"
enzsubfile <- paste("_LINCS/_KarenGuolin/", enzsubname, sep="")
enzsub <- read.table(enzsubfile, header=TRUE, sep = "\t", comment.char = "#", na.strings='', stringsAsFactors=FALSE, fill=TRUE)
# Rearrange and rename to conform to standard format
enzsub <- enzsub[,c(3,2,4,1)]
enzsub$Source <- "IP"
names(enzsub) <- c("Gene.1", "Gene.2", "Weight", "edgeType", "Source")
#_______________________________________________________________
# Note also include kinase substrate dataset from PhosphoSite
head(kinsub)
#________			
# PPI edges Functions
# This function will deal with lists containing ambiguous peptides separated by semicolons
get.gene.names.from.peps <- function(pepvec, pepsep="; ") {
     genevec=NULL
     for(i in 1:length(pepvec)) {
          x <- unlist(strsplit(as.character(pepvec[i]), pepsep, fixed=TRUE))
          genes <- unique(sapply(as.character(x),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
          genevec <- c(genevec, genes)
     }
     return(unique(genevec))
}
#
ldgenes <- get.gene.names.from.peps(rownames(ld.ratio.lim.log2))
# 4221 genes
gzgenes <- get.gene.names.from.peps(rownames(kglog2data))
# 1471 genes
# Note: Now include genes that were only detected in one experiment
gzallgenes <- unique(get.gene.names.from.peps(rownames(gzdata.all)))
# 4054 genes
length(intersect(ldgenes, gzgenes))  # 1003
length(intersect(ldgenes, gzallgenes)) # 2138
combinedgenes <- unique (c(ldgenes, gzallgenes))  # 6137
cat(combinedgenes, file=paste("_LINCS/_KarenGuolin/", "CombinedGenes.txt", sep=""), sep="\n")
# Use this to paste into Cytoscape to retrieve edges     
#....................
# useful function to find edges between nodes
filter.edges.0 <- function(nodenames, edge.file) {
     nodenames <-as.character(nodenames)
     a = as.character(edge.file$Gene.1)
     b = as.character(edge.file$Gene.2)
     edgefile.nodes <- unique(c(a,b))
     flub <- setdiff(edgefile.nodes, nodenames) 
     # show pruned nodes (turned off)
     # if (length(flub) >= 1) { 
     # cat("\n","\t", "The following GM names do not match ","\n","\t", flub) }	
     sel.edges <- edge.file[edge.file$Gene.1 %in% nodenames & edge.file$Gene.2%in% nodenames,]
     if(dim(sel.edges)[1] == 0) {return(NA)} else return(sel.edges) 
}
# Corresponding function to find all nodes connected to named nodes.
filter.edges.1 <- function(nodenames, edge.file) {
     nodenames <-as.character(nodenames)
     a = as.character(edge.file$Gene.1)
     b = as.character(edge.file$Gene.2)
     edgefile.nodes <- unique(c(a,b))
     flub <- setdiff(edgefile.nodes, nodenames) 
     # show pruned nodes (turned off)
     # if (length(flub) >= 1) { 
     #	cat("\n","\t", "The following GM names do not match ","\n","\t", flub) }
     sel.edges.1 <- edge.file[edge.file$Gene.1 %in% nodenames,]
     sel.edges.2 <- edge.file[edge.file$Gene.2%in% nodenames,]
     sel.edges <- rbind(sel.edges.1, sel.edges.2)
     if(dim(sel.edges)[1] == 0) {return(NA)} else {
          return(unique(sel.edges)) }
}
filter.ppi <- function(ppi.edges, corr.edges) {
     corr.net=corr.edges
     ppi.net=ppi.edges
     # reverse order to included reversed ppi edges
     corr.net $nodes.combined <- noquote(paste(corr.net $Gene.1, corr.net $Gene.2))
     corr.net $nodes.reversed <- noquote(paste(corr.net $Gene.2, corr.net $Gene.1))
     # make combined ppi column
     ppi.net $nodes.combined <- noquote(paste(ppi.net $Gene.1, ppi.net $Gene.2))
     # filter	ppis by corr edges	
     ppi.f1 <- ppi.net[ppi.net $nodes.combined %in% corr.net$nodes.combined, ] 
     ppi.f2 <- ppi.net[ppi.net $nodes.combined %in% corr.net$nodes.reversed, ] 
     ppi.f  <- unique(rbind(ppi.f1, ppi.f2)) 
     return(ppi.f[, 1:dim(ppi.edges)[2]])
}
# This function above can be superceded by the use of merge(); see below
# The following functions were used before the R package when STRING and GeneMANIA files were manually retrieved from the websites. But the websites have a limit on the number of nodes that can be input. So use library(STRINGdb) and the GeneMANIA Cytoscape plugin.
get.String.edgefile <- function (stringedgefile, nodenames)  {
     nodenames <- as.character(nodenames)
     stnet <- read.table(stringedgefile, header=TRUE, sep = "\t", comment.char = "", na.strings='', stringsAsFactors=FALSE, fill=TRUE)
     names(stnet)[1]="node1"  # to edit out the 'X.' or "#" 
     # fix up the gene names 
     stnet$node1 <- sapply(stnet$node1, get.gene) 
     stnet$node2 <- sapply(stnet$node2, get.gene)  
     stnodes <- unique(c(as.character(stnet$node1), as.character(stnet$node2)))
     flub <- setdiff(stnodes, nodenames) # setdiff( nodenames, stnodes) is the non-string node set
     # if there is a String flub
     if (length(flub) >= 1) { 
          cat("\n","\t", "The following String names do not match ","\n","\t", flub)
          #  Just get rid of the offending nodes
          if (any(stnet$node1 %in% flub)) stnet <- stnet[-which(stnet$node1 %in% flub), ]
          if (any(stnet$node2 %in% flub)) stnet <- stnet[-which(stnet$node2 %in% flub), ]
     }
     # 	Rearrange data file		
     stn <- stnet[,c(1:2)]
     names(stn) <- c('Gene.1', 'Gene.2')
     stne <- stn
     stne$Weight <- stnet$experimental
     stne$edgeType <- "experimental"
     stnk <- stn
     stnk$Weight <- stnet$knowledge
     stnk$edgeType <- "knowledge"
     stnh <- stn
     stnh$Weight <- stnet$homology
     stnh$edgeType <- "homology"
     stncs <- stn
     stncs$Weight <- stnet$combined_score
     stncs$edgeType <- "combined_score"
     stn <- rbind (stne, stnk, stnh, stncs)
     nones <- which (stn$Weight==0)
     stn <- stn[-nones, ]		
     return (stn)	
}
get.GM.edgefile <- function (gmfilename, nodenames) {  
     # text file from http://www.genemania.org
     nodenames <- as.character(nodenames)
     gmnet <- read.table(gmfilename, header=TRUE, sep = "\t", comment.char = "#", na.strings='', stringsAsFactors=FALSE, fill=TRUE)
     names(gmnet)[c(1:4)]=c("Gene.1", "Gene.2", "Weight", "edgeType")  
     gmnet$Gene.1 <- sapply(gmnet$Gene.1, get.gene) 
     gmnet$Gene.2 <- sapply(gmnet$Gene.2, get.gene)  
     # gmnet$Gene.1=as.factor(gmnet$Gene.1)
     # gmnet$Gene.2=as.factor(gmnet$Gene.2)
     a = as.character(gmnet$Gene.1)
     b = as.character(gmnet$Gene.2)
     gmnodes <- unique(c(a,b))
     flub <- setdiff(gmnodes, nodenames) # setdiff( nodenames, stnodes) is the non-GM node set
     # if there is an ID flub
     if (length(flub) >= 1) { 
          cat("\n","\t", "The following GM names do not match ","\n","\t", flub)
          #  Just get rid of the offending nodes
          if (any(gmnet$Gene.1 %in% flub)) gmnet <- gmnet[-which(gmnet$Gene.1 %in% flub), ]
          if (any(gmnet$Gene.2 %in% flub)) gmnet <- gmnet[-which(gmnet$Gene.2 %in% flub), ]
     }
     # 	Prune data file
     gmn <- gmnet[,c(1:4)]
     return(gmn)
}
# This function replaces some old gene names and fixes the problem that Excel converts gene names to dates (avoid Excel!).
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
#"#"#"#
#______________________________________________________________________________#
#________________________________________
#  New: library(STRINGdb)
# 		Designate local directory as /Volumes/Terra_Byte/R_Archive_2/STRINGdb
# 
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" )	
string_proteins <- string_db$get_proteins()
# dim 20457     4
test <- string_db$mp(c("FYN", "LYN", "PAG1", "CSK"))
str_sfk_Neig <- string_db$get_neighbors(test)  # 3093
string_db$plot_network(test) # cute - looks like STRING
str_sfktest <- string_db$get_interactions(test)  # 6  **
# need to convert to HUGO gene name format  
str_sfktest$Gene.1 <- string_proteins[string_proteins$protein_external_id %in% str_sfktest$from, "preferred_name"]
str_sfktest$Gene.2 <- string_proteins[string_proteins$protein_external_id %in% str_sfktest$to, "preferred_name"]
# WORKS!
# Do the whole thing
#
combined.str.ids <- string_db$mp(combinedgenes)
combined.str.net <- string_db$get_interactions(combined.str.ids)
# dim  541585     16; all: 807803     16
#	convert to HUGO gene name format 
# 	duplicates require match function 
combined.str.net$Gene.1 <- sapply(combined.str.net$from, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])
combined.str.net$Gene.2 <- sapply(combined.str.net$to, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])
# combined_score ranges from 150 to 1000
# experiments is from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID.
# 	database is from Biocarta, BioCyc, GO, KEGG, and Reactome.
#	transferred means extrapolated from other species
#    note: some of these duplcate Pathway Commons sources
dim(combined.str.net[combined.str.net$experiments>0,])
# [1] 49327    18; [2] 63871    18
dim(combined.str.net[combined.str.net$experiments_transferred>0,])
# [1] 195352    18; [2] 298450     18
str.e <- combined.str.net[combined.str.net$experiments>0,]
str.et <- combined.str.net[combined.str.net$experiments_transferred>0,]
str.d <- combined.str.net[combined.str.net$database>0,]
str.dt <- combined.str.net[combined.str.net$database_transferred>0,]
combined.str <- unique(rbind (str.e, str.et, str.d, str.dt))  # 251613 edges; 374352 v2
combined.str$edgeType <- "STRINGdb"
combined.str[combined.str$database>0, "edgeType"] <- "database"
combined.str[combined.str$database_transferred>0, "edgeType"] <- "database"
combined.str[combined.str$experiments>0, "edgeType"] <- "experiments"
combined.str[combined.str$experiments_transferred>0, "edgeType"] <- "experiments"
combined.str$Weight <- rowSums(combined.str[, c("experiments", "experiments_transferred", "database", "database_transferred")])
range(combined.str$Weight)  # 43 to 2859
range(log(combined.str$Weight))  #  3.761200 7.958227
combined.str.edges <- combined.str[,c("Gene.1", "Gene.2", "Weight", "edgeType")]
# To make them on the same scale as other dbs:
combined.str.edges$Weight <- combined.str.edges$Weight/1000	
# ****  edges range(combined.str.edges$Weight) 0.043 2.859
save(ldgenes, gzallgenes, combinedgenes, combined.str, combined.str.edges, combined.str.ids, combined.str.net,  file=paste(pcpath, "GZ_PPI_Networks2.RData", sep=""))
#####-------------------------	
#  GeneMANIA Cytoscape app has the ability to export network as text in the Results panel. Manually duplcate the file and delete all but the PPIs for
gmfilename <- "Combined_GeneMANIA_PPIs2.txt"
#cat(combinedgenes)
# 
combined.GM.net <- read.table(paste(pcpath, gmfilename, sep=""), header=TRUE, sep = "\t", comment.char = "#", na.strings='', quote = "", stringsAsFactors=FALSE, fill=TRUE)
# dim 567994  6; [2] 1009527      5
# combined.GM.net <- combined.GM.net[, 1:5]
# combined.GM.edges <- get.GM.edgefile(paste(pcpath, gmfilename, sep=""), combinedgenes)
# dim  550380      4; [2] 686055      4
#  	 The following GM names do not match  
# CCNA2 CCNB1 CDK4 CKS2 ILK RPA3 RPL32 SRPK1 TUBG1 UPF2 VCAM1 CUL4A ETF1 NHP2 PJA2 POLR2E rcl SNRPD2 SUZ12 UBC RNF2
# range(combined.GM.net$Weight)
# [1] 7.509860e-08 6.920997e-01
# last edge ZRANB2  TRA2B
# Note: in Cytoscape only the following edges should be included; otherwise edit manually or here.
unique(combined.GM.net$Type)
#  "Genetic Interactions"  "Pathway" "Physical Interactions" "Predicted"   
# Note also that Genetic Interactions may be indirect, and may be excluded for some purposes.
# cbGM <- combined.GM.edges[-grep("Shared protein domains", combined.GM.edges$edgeType)]
combined.GM.edges <- combined.GM.net
names(combined.GM.edges)[4] <- "edgeType"
# *****
# Now do Pathway Commons and Bioplex	
combined.pc.edges <- filter.edges.0(combinedgenes, pcnet)  # 119255 edges	
combined.bp.edges <- filter.edges.0(combinedgenes, bioplex)  # 6328 edges
# combined.bp.edges <- combined.bp.edges[,c(1,2,4,3)]
combined.enzsub.edges <- filter.edges.0(combinedgenes, enzsub)  # 4210 edges	
# 
combined.kinsub.edges <- filter.edges.0(combinedgenes, kinsub) # 4086 edges; 4428 with 2019 kinsub
# Make and "all ppi edge" file for comparisone
combined.all.ppi <- rbind(combined.enzsub.edges[,1:4], combined.pc.edges[,1:4], combined.bp.edges, combined.GM.edges[,1:4], combined.str.edges, combined.kinsub.edges)
# 1572925 edges!  1577353 with new kinsub
combined.all.ppi.g <- graph_from_data_frame(combined.all.ppi)
#   927576 edges; 1009527 with kinsub
#-------------------------------------------------------------------------
# NEXT: Cluster-filtered networks
# 	We want PPI edges, filtered to show only interactions between proteins whose modifications co-clustered

# 		Alternative method to filter.ppi() function
unique(combined.all.ppi$edgeType)
head(combined.all.ppi[combined.all.ppi$edgeType=="RnaReference",]) # only two - auto loops
combined.all.ppi <- combined.all.ppi[-which(combined.all.ppi$edgeType=="RnaReference"),]
head(combined.all.ppi[combined.all.ppi$edgeType=="catalysis-precedes",]) # 14648 leave in for now
# Note that some of these are directed
diredgeTypes <- c("pp", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of",  "controls-state-change-of", "ACETYLATION", "METHYLATION", "PHOSPHORYLATION", "SUMOYLATION", "UBIQUITINATION", "catalysis-precedes") 
gzalltgenecccn.edges2 <- gzalltgenecccn.edges[, c("Gene.2", "Gene.1", "Weight", "interaction")]
names(gzalltgenecccn.edges2) <- c("Gene.1", "Gene.2", "Weight", "interaction") 
# Correlation goes both ways. 
gzalltgenecccn.edges3 <- rbind(gzalltgenecccn.edges, gzalltgenecccn.edges2)
# need only Gene.Names
gzalltgene.cfn.genes <- gzalltgenecccn.edges3[,1:2]
##>>>>
# Now Filter 
gzalltgene.cfn <- merge(gzalltgene.cfn.genes, combined.all.ppi, by=c("Gene.1", "Gene.2"))
# dim 13591     13664 with kinsub; 13671 with 2019 kinsub
gzalltgene.cfn.g <- graph_from_data_frame(gzalltgene.cfn)
# save(ldgenes, gzallgenes, combinedgenes, combined.str, combined.str.edges, combined.str.ids, combined.str.net, combined.GM.net, combined.GM.edges, combined.pc.edges, combined.bp.edges, combined.enzsub.edges, combined.all.ppi, combined.all.ppi.g, gzalltgene.cfn, gzalltgene.cfn.g, file=paste("_LINCS/_KarenGuolin/", "GZ_PPI_Networks2.RData", sep=""))
#-------------------------------------------------------------------------
#---------+++===>>> Plot betweenness to compare different cfns.
# Start with same combined.all.ppi for ld data
# 
load(file="/Users/Mark_Grimes/Dropbox/_Work/R_/_LINCS/Bin_Curated_Data/ld_tsne5.RData") 
# Contains essl.gene.cccn0; or: 
essl.gene.cccn0 <- read.table("/Users/Mark_Grimes/Dropbox/_Work/R_/_LINCS/Bin_Curated_Data/lim_log2_clusters_Genes_CCCN0.txt", header=TRUE, sep = "\t", comment.char = "#", na.strings='', stringsAsFactors=FALSE, fill=TRUE)
# Make igraph object
ld.gene.cccn.g <- graph.adjacency(as.matrix(essl.gene.cccn0), mode="lower", diag=FALSE, weighted="Weight")
#  
hist(edge_attr(ld.gene.cccn.g)[[1]], breaks=1000, col="deepskyblue", ylim=c(0,1200))
# to get the adjacency matrix back, but NOTE, edges are 0 or 1.
# gzallt.gene.cccn.mat <- as_adjacency_matrix(gzallt.gene.cccn.g, type="both")
# ___make an edge list file
ldgenegenecccn.edges <- data.frame(as_edgelist(ld.gene.cccn.g))
names(ldgenegenecccn.edges) <- c("Gene.1", "Gene.2")
ldgenegenecccn.edges$Weight <- edge_attr(ld.gene.cccn.g)[[1]]
ldgenegenecccn.edges$interaction <- "correlation" 
ldgenegenecccn.edges $interaction[ldgenegenecccn.edges$Weight<=-0.5] <- "negative correlation"
ldgenegenecccn.edges $interaction[ldgenegenecccn.edges$Weight>=0.5] <- "positive correlation"
# dim 35243      4
ldgenegenecccn.edges2 <- ldgenegenecccn.edges[, c("Gene.2", "Gene.1", "Weight", "interaction")]
names(ldgenegenecccn.edges2) <- c("Gene.1", "Gene.2", "Weight", "interaction") 
# Correlation goes both ways. 
ldgenegenecccn.edges3 <- rbind(ldgenegenecccn.edges, ldgenegenecccn.edges2)
# need only Gene.Names
ldgene.cfn.genes <- ldgenegenecccn.edges3[,1:2]
# Now Filter 
ldgene.cfn <- merge(ldgene.cfn.genes, combined.all.ppi, by=c("Gene.1", "Gene.2"))
# dim 7967     8000 with kinsub; 2 more with 2019 kinsub
ldgene.cfn.g <-  graph.data.frame(ldgene.cfn, directed=FALSE)



#--------------------
# Explore network attributes
#
make.netatts <- function(ig.network, ig.unfiltered, keyfile, ppinet=FALSE){
        cccndegree  <- igraph::degree(ig.network, mode="all", loops=F, normalized=F)
        cccnbetween <- igraph::betweenness(ig.network)
        net.df <- data.frame(cccndegree, cccnbetween)            
        net.df$Peptide.Name <- rownames(net.df) 
        net.df <- net.df[,c(3,1,2)]
        if(ppinet==TRUE){
                names(net.df) <- c("Gene.Name", "ppidegree", "ppibetween")
                allppibetween <- betweenness(ig.unfiltered)
                allppibetween <- allppibetween[sort(names(allppibetween))]
                allbetween.df <- data.frame(allppibetween)
                allbetween.df$Gene.Name <- rownames(allbetween.df)
                net.df <- merge(net.df, allbetween.df, all=TRUE)
                net.df$norm.ppibetween <- net.df$ppibetween/net.df$allppibetween}
        netatts.df <- merge(keyfile, net.df, all=TRUE)
        netatts.df[is.na(netatts.df)] <- 0
        return(netatts.df)
}

gzgene.cfn.netatts <- make.netatts(ig.network=gzalltgene.cfn.g, ig.unfiltered=combined.all.ppi.g, keyfile=gzallt.gene.key, ppinet=TRUE)
# ld.gene.key is in "ld_networks3.RData"
load(file="/Users/Mark_Grimes/Dropbox/_Work/R_/_LINCS/Bin_Curated_Data/ld_networks3.RData") 
ldgene.cfn.netatts <- make.netatts(ig.network=ldgene.cfn.g, ig.unfiltered=combined.all.ppi.g, keyfile=ld.gene.key, ppinet=TRUE)

compare.nettatts <- merge(gzgene.cfn.netatts, ldgene.cfn.netatts, by="Gene.Name", all=F)
# 6167 out of 6188/6187

plot(compare.nettatts$ppibetween.y~compare.nettatts$ppibetween.x, col="blue", pch=17)
summary(lm(compare.nettatts$ppibetween.y~compare.nettatts$ppibetween.x))
# R-squared:  0.271
plot(compare.nettatts$ppidegree.y~compare.nettatts$ppidegree.x, col="red", pch=16)
summary(lm(compare.nettatts$ppidegree.y~compare.nettatts$ppidegree.x))
# R-squared:  0.3673
# prune to those with data in both
com.natts <- compare.nettatts[which(compare.nettatts$No.Samples.x>0 & compare.nettatts$No.Samples.y>0),]
# now 1878 
plot(com.natts$ppibetween.y~com.natts$ppibetween.x, col="blue", pch=17)
summary(lm(com.natts$ppibetween.y~com.natts$ppibetween.x))
# R-squared:  0.2786
plot(com.natts$ppidegree.y~com.natts$ppidegree.x, col="red", pch=16)
summary(lm(com.natts$ppidegree.y~com.natts$ppidegree.x))
 # R-squared:  0.3662

write.table(ldgene.cfn.netatts, file=paste("_LINCS/_KarenGuolin/","revised_ldgene_attributes.txt", sep=""), sep="\t", eol = "\n", quote = FALSE, row.names=TRUE, col.names=TRUE)
write.table(gzgene.cfn.netatts, file= paste("_LINCS/_KarenGuolin/","GZ_allt_cfn_attributes.txt", sep=""), sep="\t", eol = "\n", quote = FALSE, row.names=TRUE, col.names=TRUE)
##############################################################################################
# Examine subsets of edges
edgeslistname <- "/Users/Mark_Grimes/Dropbox/_Work/R_/_LINCS/_KarenGuolin/ten_cell_lines_cfn_edge_types.txt"
edgeslist <- read.table(edgeslistname, sep = "\t", skip = 0, header=TRUE, blank.lines.skip=T, fill=T, quote="\"", dec=".", comment.char = "", stringsAsFactors=F)
gzedgeTypes <- unique(gzalltgene.cfn$edgeType)
edgeslist[grep('physical', edgeslist$Include.in), "CFN.Edge.Type"]
physicaledges <- c("ACETYLATION", "METHYLATION", "PHOSPHORYLATION", "pp", "experiments", "in-complex-with", "interacts-with", "Physical Interactions", "BioPlex")
gzalltgene.physical.cfn <- gzalltgene.cfn[gzalltgene.cfn$edgeType %in% physicaledges,]
# 8398/13591 edges (or 8471/13664 with kinsub; 7 more 2019 kinsub)
ldedgeTypes <- unique(ldgene.cfn$edgeType)
ldgene.physical.cfn <- ldgene.cfn[ldgene.cfn$edgeType %in% physicaledges,]
# 4912/7967 edges;  4945/8000 with kinsub; 2 more 2019 kinsub
# What is the overlap?
# Convert to new naming method
# names(gzalltgene.physical.cfn) <- names(ldgene.physical.cfn)
overlap.cfn <- merge(ldgene.physical.cfn, gzalltgene.physical.cfn, by=c("Gene.1", "Gene.2"), all=FALSE)
# 1201 edges; 1261 with kinsub + 1 more 2019 kinsub
# NOTE: this contains duplicated edges, see:
filter.edges.0(c("ACTB","CFL1"), overlap.cfn)
# Test if some are in reverse orientation (binding interactions only, not enzymatic)
#>>>>>
# NOTE: Reversing nodes manually also duplicates edges
        ld.p.cfn.rev <- ldgene.physical.cfn[ldgene.physical.cfn$edgeType %in% c("experiments", "in-complex-with", "interacts-with", "Physical Interactions", "BioPlex"), c(2,1,3,4)]
        names(ld.p.cfn.rev) <- names(ldgene.cfn)
        gz.p.cfn.rev <- gzalltgene.physical.cfn[gzalltgene.physical.cfn$edgeType %in% c("experiments", "in-complex-with", "interacts-with", "Physical Interactions", "BioPlex"), c(2,1,3,4)]
        names(gz.p.cfn.rev) <- names(ldgene.cfn)
        overlap.cfn2 <- merge(rbind(ldgene.physical.cfn, ld.p.cfn.rev), rbind(gzalltgene.physical.cfn, gz.p.cfn.rev), by=c("Gene.1", "Gene.2"), all=FALSE)
        overlap.cfn2$Weight <- rowSums(overlap.cfn2[, grep("Weight", names(overlap.cfn2))])
        overlap.cfn2$edgeType <- sapply(overlap.cfn2[, grep("edgeType", names(overlap.cfn2))], function(x) paste(c(x), collapse=", "))
        #apply(overlap.cfn2[, grep("edgeType", names(overlap.cfn2))], 1, function(x) paste(unlist(unique(x))))
        overlap.cfn2 <- overlap.cfn2[,c(1,2, 7,8)]
#>>>>
 # New Approach: Sort the node names to ensure undirected edges are the same 
        physicaledges <- c("ACETYLATION", "METHYLATION", "PHOSPHORYLATION", "pp", "experiments", "in-complex-with", "interacts-with", "Physical Interactions", "BioPlex")
        unique(gzalltgene.physical.cfn$edgeType)
        unique(ldgene.physical.cfn$edgeType) # no ACETYLATION
        # Separate directed and undirected edges        
        ld.undir.sub <- ldgene.physical.cfn[ldgene.physical.cfn$edgeType %in% physicaledges[5:9], ] #4872
        gz.undir.sub <- gzalltgene.physical.cfn[gzalltgene.physical.cfn$edgeType %in% physicaledges[5:9],]# 8314
        ld.dir.sub <- ldgene.physical.cfn[ldgene.physical.cfn$edgeType %in% physicaledges[1:4], ] # 75
        gz.dir.sub <- gzalltgene.physical.cfn[gzalltgene.physical.cfn$edgeType %in% physicaledges[1:4],] # 164
# Sort node names - for undirected edges only
sort.nodes <- function(edgefile) {
        sorted.df <- t(apply(edgefile[, 1:2], 1, sort))
        newfile <- cbind(sorted.df, edgefile[, 3:ncol(edgefile)])
        names(newfile)[1:2] <- names(edgefile)[1:2]
        return(newfile)
}
        edgefile <- ld.undir.sub
        test <- rbind(head(edgefile), tail(edgefile))
look.2 <- sort.nodes(test)        
 # Okay, works
ld.undir.sub.sort <- sort.nodes(ld.undir.sub)        
gz.undir.sub.sort <- sort.nodes(gz.undir.sub)   
overlap.cfn2 <- merge(rbind(ld.dir.sub, ld.undir.sub.sort), rbind(gz.dir.sub, gz.undir.sub.sort), all=FALSE)
# NOTE 1558 edges with by=c("Gene.1", "Gene.2"), some duplicated; only 508 with by = intersect(names(x), names(y))
#________________________________________________
# igraph objects
ldgene.physical.cfn.g <-  graph.data.frame(ldgene.physical.cfn, directed=FALSE)
gzalltgene.physical.cfn.g <-  graph.data.frame(gzalltgene.physical.cfn, directed=FALSE)
overlap.cfn2.g <-  graph.data.frame(overlap.cfn2, directed=FALSE)


#--------------------------------------------------------------------------------------------------
# Make a simplifed version of the CFN with edges merged using mergeEdges.RCy32 (see New_RCy3_Functions.R)

# Make edgefiles RCy3 format with edgeType.to.interaction 
#--Make new edgefiles
ldgene.cfn.rcy3 <- edgeType.to.interaction(ldgene.cfn)
ldgene.physical.cfn.rcy3 <- edgeType.to.interaction(ldgene.physical.cfn)
gzalltgene.cfn.rcy3 <- edgeType.to.interaction(gzalltgene.cfn)
gzalltgene.physical.cfn.rcy3 <- edgeType.to.interaction(gzalltgene.physical.cfn)
overlap.cfn2.rcy3 <- edgeType.to.interaction(overlap.cfn2)
ldgene.physical.cfn.merged <- mergeEdges.RCy32(ldgene.physical.cfn.rcy3)
gzalltgene.physical.cfn.merged <- mergeEdges.RCy32(gzalltgene.physical.cfn.rcy3)
overlap.cfn2.merged <- mergeEdges.RCy32(overlap.cfn2.rcy3) # now 254 edges
# igraph objects
ldgene.physical.cfn.merged.g <-  graph.data.frame(ldgene.physical.cfn.merged, directed=FALSE)
gzalltgene.physical.cfn.merged.g <-  graph.data.frame(gzalltgene.physical.cfn.merged, directed=FALSE)
overlap.cfn2.merged.g <-  graph.data.frame(overlap.cfn2.merged, directed=FALSE)
#__________________________________________
# From KGnegcorredges2.R
gzallt.cccnplus <- gzallt.cccn.edges.plus[,1:4]
names(gzallt.cccnplus) <- c("source", "target", "Weight", "interaction")
gzallt.gpe <- gzalltgenepep.edges[,c(1,2,4,3)]
names(gzallt.gpe) <- c("source", "target", "Weight", "interaction")
gzallt.gpe$Weight <- as.numeric(as.character(gz.gpe$Weight))
# Prune gpe nodes not in cccn
cccn.nodes <- unique(c(gzallt.cccnplus$source, gzallt.cccnplus$target)) # 7715
cfn.nodes <- unique(c(gz.cfn$source, gz.cfn$target))
gpe.nodes <- unique(gz.gpe$source) # 1986
gzallt.gpe.pruned <- gzallt.gpe[gzallt.gpe$source %in% cccn.nodes,]
gzallt.gpe.pruned <- gzallt.gpe.pruned[gzallt.gpe.pruned$target %in% cfn.nodes,]
gzallt.network <- rbind(gz.cfn, gzallt.gpe.pruned, gzallt.cccnplus)
dim(gzallt.network)  # [1] 95160     4
## __________________________________________
# Incorportate physcial CFNs from PPI_Networks.R
# Edgefile = gzalltgene.physical.cfn.merged
gzallt.physical.network <- rbind(gzalltgene.physical.cfn.merged, gzallt.gpe.pruned, gzallt.cccnplus)
# 86421 edges
gzallt.physical.network <- edgeType.to.interaction(gzallt.physical.network)

save(ldgenes, gzallgenes, combinedgenes, combined.str, combined.str.edges, combined.str.ids, combined.str.net, combined.GM.net, combined.GM.edges, combined.pc.edges, combined.bp.edges, combined.enzsub.edges, combined.all.ppi, combined.all.ppi.g, gzallt.key, gzallt.cccnplus, gzalltgene.cfn, gzalltgene.cfn.g, gzallt.gene.key, gzallt.network, gzallt.network.all, gzallt.physical.network, gzalltgene.physical.cfn.merged, gzalltgene.physical.cfn.merged.g, gzallt.cccnplus, gzalltgene.ave.ratios, gz.cfn, gzallt.network.all, gzallt.all.cf, gzallt.all.cf.pruned, gzalltgene.cf, gzalltgene.all.cf, gzalltgene.all.cf.pruned, gz.cf, gz.cf.pruned, ld.gene.key, ld.gene.cccn.g, ldgene.cfn, ldgene.cfn.g, gzgene.cfn.netatts, ldgene.cfn.netatts, compare.nettatts, ldgene.physical.cfn, gzalltgene.physical.cfn, overlap.cfn, overlap.cfn2, ldgene.physical.cfn.g, gzalltgene.physical.cfn.g, overlap.cfn2.g, gzalltgene.cfn.rcy3, gzalltgene.physical.cfn.rcy3, ldgene.cfn.rcy3, ldgene.physical.cfn.rcy3, overlap.cfn2.rcy3, ldgene.physical.cfn.merged, ldgene.physical.cfn.merged.g, overlap.cfn2.merged, overlap.cfn2.merged.g,  keggpath, wikipath, bioplanet, file=paste("_LINCS/_KarenGuolin/", "GZ_PPI_Networks2.RData", sep=""))





#-----------------------------------------------------------------------------------------------
ldgene.cfn.cf <- ldgene.cfn.netatts
ldgene.cfn.cf$id <- ldgene.cfn.cf$Gene.Name
ldgene.cfn.cf <- ldgene.cfn.cf[, c(8,1:7)]
names(ldgene.cfn) <- c("source", "target", "Weight", "interaction")

ldnetwork.suid <- createNetworkFromDataFrames(ldgene.cfn.netatts, ldgene.cfn, title="Revised LD CFN", collection = "Interactions")
setLayoutProperties("genemania-force-directed", list(numIterations=100, defaultSpringCoefficient=0.1, defaultSpringLength=50, minNodeMass=0.001, maxNodeMass=2000, midpointEdges=250, curveSteepness=7.0e-03, isDeterministic=1, singlePartition=0, ignoreHiddenElements=1))
layoutNetwork("genemania-force-directed")


#-----------------------------------------------------------------------------------------------
# Cuneyt's Network from cluster150perc0dot7.txt
nodedata <- getTableColumns("node", columns = c("name", "Column 4"))
unique(nodedata$`Column 4`)
# [1] "Gene" ""     "P"   
edgedata <- getTableColumns("edge")
names(edgedata)
# This is a series of Genes interacting with Pathways, e.g.,
# [997] "PI3K-Akt Signaling Pathway WP4172 (interacts with) PTEN"    
# Turn it into an edgefile
edgenames <- edgedata$name
# Switch so Gene is source
target <- as.character(sapply(edgenames, function (x) unlist(strsplit(x, " (", fixed=TRUE), use.names=F)[1]))
source <- as.character(sapply(edgenames,  function (x) unlist(strsplit(x, ") ", fixed=TRUE), use.names=F)[2]))
pathedges <- data.frame(source, target, Weight=1, interaction="in pathway")
filter.edges.1.RCy3   ("MAPK1", pathedges)      # 23 pathways  
filter.edges.1.RCy3   ("ROS1", pathedges)       # only 1 
filter.edges.1.RCy3   ("HSP90AA1", pathedges)        # 4
filter.edges.1.RCy3   ("EP300", pathedges)        # 6
# Can we look the other way around?
WntPathGenes <- filter.edges.1.RCy3 ("Wnt Signaling Pathway and Pluripotency WP399", pathedges)$source
WntGenes <- WntPathGenes[WntPathGenes %in% c(ldgenes, gzallgenes)]  # 46/104 **
WntGeneGZNet <- filter.edges.1.RCy3(WntGenes, gzalltgene.physical.cfn.merged)  # 215     
WntGeneLDNet <- filter.edges.1.RCy3(WntGenes, ldgene.physical.cfn.merged)       # 53   
WntGeneOverlapNet <- filter.edges.1.RCy3(WntGenes, overlap.cfn2.merged)         # 5 edges
WntOverlapGenes <- extract.gene.names.RCy3(WntGeneOverlapNet)
wntoverlap.cf <-  gz.cf[gz.cf$Gene.Name %in% WntOverlapGenes,]
ros1ep300 <- connectNodes.all.RCy3.exclude(c("ROS1", "EP300"), ig.graph=gzalltgene.physical.cfn.merged.g, edgefile=gzalltgene.physical.cfn.merged, exclude = "MYH9")
# ROS1 is not well represented enough in the data to be part of any clusters
        filter.edges.1.RCy3("ROS1", gzalltgene.physical.cfn.merged) # NA
        filter.edges.1.RCy3("ROS1", ldgene.physical.cfn.merged) # NA
        filter.edges.1("ROS1", essl.cfn) # NA
        ros1ldcccn <- essl.cccn.edges[grep("ROS1", essl.cccn.edges$Gene.1) | grep("ROS1", essl.cccn.edges$Gene.2),]
        ros1gzcccn <- gzgenecccn.edges[grep("ROS1", gzgenecccn.edges$Gene.1) | grep("ROS1", gzgenecccn.edges$Gene.2),]
        #NADA
# BUT: ROS1 is in the combined CFN/CCCN:
        filter.edges.1.RCy3("ROS1", gzallt.physical.network) 
        # ZNF24
        filter.edges.1.RCy3("ZNF24", gzallt.physical.network) 
        
# Pathways
        # 
 #       KEGG_2019_Human.txt, WikiPathways_2019_Human.txt
        # New: bioplanet_pathway.csv
bioplanetname <- paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/",  "bioplanet_pathway.csv", sep="")
bioplanet <- read.csv(bioplanetname, stringsAsFactors = F)
dim(bioplanet)
# 74148     4
length(unique(bioplanet$PATHWAY_ID))
# 1658
# Separtate into a list as for kegg and wikipath
bp <- dlply(bioplanet, .(PATHWAY_NAME))
# Remove the first two vector elements from each list element (iteratively)
bpy <- lapply(bp, `[`, -c(1:3))
bpz <- lapply(bpy, unlist, use.names = FALSE)
#y <- lapply(y, function(x) x[-1]) # same as above for first one
bioplanet <- bpz 
# KEGG pathways
keggpath <- y
keggname <- paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/",  "KEGG_2019_Human.txt", sep="")
wikiname <- paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/",  "WikiPathways_2019_Human.txt", sep="")
keggpath <- read.table(keggname, header=TRUE, sep = "\t", comment.char = "#", na.strings='', stringsAsFactors=FALSE, fill=TRUE)#  
# Read table doesn work try another methods to # Read in the data
keggpath <- scan(keggname, what="", sep="\n")
# Separate elements by one or more whitepace
# y <- strsplit(keggpath, "[[:space:]]+")
z <- strsplit(keggpath, "\t")
# Extract the first vector element and set it as the list element name
names(z) <- sapply(z, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first two vector elements from each list element (iteratively)
y <- lapply(z, `[`, -c(1:2))
#y <- lapply(y, function(x) x[-1]) # same as above for first one
keggpath <- y
 # Wikipathways
wikipath <- scan(wikiname, what="", sep="\n")
w <- strsplit(wikipath, "\t")
# Extract the first vector element and set it as the list element name
names(w) <- sapply(w, `[[`, 1)
wikipath <- lapply(w, `[`, -c(1:2))
# See who is where
keggpath[grep("ROS1", keggpath)]
wikipath[grep("ROS1", wikipath)]
#:# What is pathway intersection
a <- wikipath[["Spinal Cord Injury WP2431"]]
b <- wikipath[["VEGFA-VEGFR2 Signaling Pathway WP3888"]]
ab <- intersect (a, b)
abgz <- ab[ab %in% gzallt.gene.key$Gene.Name]
filter.edges.0.RCy3(abgz, gzalltgene.physical.cfn) # NA ??
filter.edges.0.RCy3(abgz, ldgene.physical.cfn) # NA ??
test <- composite.shortest.paths(genes1=abgz, genes2=abgz, network=gzalltgene.physical.cfn.merged, exclude="MYH9")
# 468 edges; 204 with merged
test.2 <- composite.shortest.paths(genes1=a, genes2=b, network=gzalltgene.physical.cfn.merged, exclude="MYH9")


#_________________________________________
```{r plot simpler sub-network}
# The simplest way to get the node attributes is to select the genes in the entire network above and the node attribute file.
genenames <- extract.gene.names.RCy3(alkep300)
alkep300.cf <- gz.cf[gz.cf$Gene.Name %in% genenames,]
# Note: this contains PTMs. If you want only the CFN (gene nodes) use:
alkep300.gene.cf <- alkep300.cf[which(alkep300.cf$Node.ID=="gene"),]
