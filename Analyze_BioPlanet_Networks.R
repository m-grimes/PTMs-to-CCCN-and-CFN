# Analysis of BioPlanet Pathway Interaction Networks
# July 9, 2021
# Mark Grimes
#
# Networks created in in "Create_BioPlanet_Network.R"
#============================================================================== #
localpath = getwd()	
comp_path <- unlist(strsplit(localpath, "Dropbox"))[[1]]
if(length(grep("Dropbox", comp_path))==0) { 
  comp_path <- unlist(strsplit(localpath, "Documents"))[[1]]}
source(paste(comp_path, "/Dropbox/_Work/R_/MG_packages.R", sep=""))
# load work
load(file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "GZ_PPI_Networks2.RData", sep=""))
load(file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "TenCell.RData", sep=""))
load(file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "BioPlanetNetworks.RData", sep=""))
# Source drug groupings script
source(file=paste(comp_path, "/Documents/PTMs-to-CCCN-and-CFN/", "Drug Groupings.R", sep=""))
#
# Pathway-pathway networks:
# total.pathway.net: Weight = Weight.clust - Weight.bp; Combined.Weight = Weight.clust + Weight.bp
# total.pathway.net.no.bp: Weight.bp = 0; Weight = Weight.clust
# pathway.crosstalk.network: Weight = either cluster evidence or Jaccard similarity (from bioplanet) as individual edges
#______________________________________________________________________________
# What is the size of the network without any threshold?
#  dim (total.pathway.net) 645709      4
# What does the entire pathway crostalk network graph look like? 
no.nodes.total.pathway.net <- length(unique(c(total.pathway.net$source, total.pathway.net$target))) # 1488
# no. possible edges = 0.5*no.nodes*(no.nodes-1)
no.poss.total.edges <- 0.5*(no.nodes.total.pathway.net*(no.nodes.total.pathway.net-1))
# 1106328
total.pathway.net.density <- dim (total.pathway.net)[1]/no.poss.total.edges
# density is now 0.5836506
#-------------------------------------------------------------------------------
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
# What does the entire pathway crostalk network graph look like? Which edge weights to focus on? 
plot(total.pathway.net$Weight.bp~total.pathway.net$Weight.clust)
plot(total.pathway.net$Combined.Weight~total.pathway.net$Weight.clust)
plot(total.pathway.net$Weight.normalized~total.pathway.net$Weight.clust)
dev.new()
plot(total.pathway.net$Weight.normalized~total.pathway.net$Weight.clust, ylab= "Transformed weight", xlab="Cluster evidence weight", pch=20, col=alpha(col2hex("green"), 0.1), log="xy", cex=1.8)
points(total.pathway.net$Combined.Weight~total.pathway.net$Weight.clust,  pch=20, col=alpha(col2hex("magenta"), 0.1), cex=1.8)
points(big.path.net$Weight~big.path.net$Weight.clust,  pch=20, col=alpha(col2hex("blue"), 0.1), cex=1.8)
# *** ClusterBPWeightGraph2.pdf
legend("bottomright", pt.cex=1.8, pch=20, col=c("green", "magenta", "blue"), legend = c("Weight normalized", "Weight combined", "Weight normalized, filtered"))
# Without log
dev.new()
plot(total.pathway.net$Weight.normalized~total.pathway.net$Weight.clust, ylab= "Transformed weight", xlab="Cluster evidence weight", pch=20, col=alpha(col2hex("green"), 0.25), cex=1.8)
points(total.pathway.net$Combined.Weight~total.pathway.net$Weight.clust,  pch=20, col=alpha(col2hex("magenta"), 0.2), cex=1.8)
points(big.path.net$Weight~big.path.net$Weight.clust,  pch=20, col=alpha(col2hex("blue"), 0.2), cex=1.8)
# *** ClusterBPWeightGraph2lin.pdf
legend("bottomright", pt.cex=1.8, pch=20, col=c("green", "magenta", "blue"), legend = c("Weight normalized", "Weight combined", "Weight normalized, filtered"))
# *** ClusterBPWeightGraph2lin.pdf
# Without log: normalized
dev.new()
plot(tpnn$Weight.normalized~tpnn$Weight.clust, ylab= "Transformed weight", xlab="Cluster evidence weight", pch=20, col=alpha(col2hex("green"), 0.25), cex=1.8)
points(tpnn$Combined.Weight~tpnn$Weight.clust,  pch=20, col=alpha(col2hex("magenta"), 0.2), cex=1.8)
points(tpnn[which(tpnn$Weight.bp==0), "Combined.Weight"]~tpnn[which(tpnn$Weight.bp==0), "Weight.clust"],  pch=20, col=alpha(col2hex("blue"), 0.5), cex=1.8)
legend("bottomright", pt.cex=1.8, pch=20, col=c("green", "magenta", "blue"), legend = c("Weight normalized", "Weight combined", "No common genes"))
# *** ClusterBPWeightGraph2lin3.pdf
# Normalized, with log
dev.new()
plot(tpnn$Weight.normalized~tpnn$Weight.clust, ylab= "Transformed weight", xlab="Cluster evidence weight", pch=20, col=alpha(col2hex("green"), 0.25), log="xy", cex=1.8)
# Warning: 162514 y values <= 0 omitted from logarithmic plot
points(tpnn$Combined.Weight~tpnn$Weight.clust,  pch=20, col=alpha(col2hex("magenta"), 0.2), cex=1.8)
points(tpnn[which(tpnn$Weight.bp==0), "Combined.Weight"]~tpnn[which(tpnn$Weight.bp==0), "Weight.clust"],  pch=20, col=alpha(col2hex("blue"), 0.2), cex=1.8)
legend("bottomright", pt.cex=1.8, pch=20, col=c("green", "magenta", "blue"), legend = c("Weight normalized", "Weight combined", "No common genes"))
# *** ClusterBPWeightGraph2log3.pdf
dev.new()
plot(total.pathway.net$Weight.bp~total.pathway.net$Weight.clust, ylab= "Bioplanet Jaccard similarity", xlab="Cluster evidence weight", pch=20, col=alpha(col2hex("green4"), 0.25), cex=1.8)
summary(lm(total.pathway.net$Weight.bp~total.pathway.net$Weight.clust))
# R-squared:  0.02037 
points(total.pathway.net[which(total.pathway.net$Weight.clust==0), "Weight.bp"]~total.pathway.net.no.bp[which(total.pathway.net$Weight.clust==0), "Weight.clust"],  pch=20, col=alpha(col2hex("red"), 0.2), cex=1.8) 
#none, duh
points(total.pathway.net[which(total.pathway.net$Weight.bp<=max(total.pathway.net$Weight.bp)/10), "Weight.bp"]~total.pathway.net.no.bp[which(total.pathway.net$Weight.bp<=max(total.pathway.net$Weight.bp)/10), "Weight.clust"],  pch=20, col=alpha(col2hex("red"), 0.25), cex=1.8)
points(total.pathway.net[which(total.pathway.net$Weight.bp==0), "Weight.bp"]~total.pathway.net.no.bp[which(total.pathway.net$Weight.bp==0), "Weight.clust"],  pch=20, col=alpha(col2hex("darkblue"), 0.2), cex=1.8)
legend("topright", pt.cex=1.8, pch=20, col=c("green4", "red", "darkblue"), legend = c("Weights compared", "Weight bioplanet < 10%", "Weight bioplanet = 0"))
# *** Weight.clust.vs.weight.bp.pdf/png
# Normalizde version
dev.new()
plot(tpnn$Weight.bp~tpnn$Weight.clust, ylab= "Bioplanet Jaccard similarity", xlab="Cluster evidence weight", pch=20, col=alpha(col2hex("green4"), 0.25), cex=1.8)
summary(lm(tpnn$Weight.bp~tpnn$Weight.clust))
# R-squared:  0.02037 (unchanged)
points(tpnn[which(tpnn$Weight.bp<=max(tpnn$Weight.bp)/10), "Weight.bp"]~tpnn[which(tpnn$Weight.bp<=max(tpnn$Weight.bp)/10), "Weight.clust"],  pch=20, col=alpha(col2hex("red"), 0.25), cex=1.8)
points(tpnn[which(tpnn$Weight.bp==0), "Weight.bp"]~tpnn[which(tpnn$Weight.bp==0), "Weight.clust"],  pch=20, col=alpha(col2hex("darkblue"), 0.25), cex=1.8)
legend("topright", pt.cex=1.8, pch=20, col=c("green4", "red", "darkblue"), legend = c("Weights compared", "Weight bioplanet < 10%", "Weight bioplanet = 0"))
# Weight.clust.vs.weight.bp2.pdf
# Conclusions: 
# 1. There is little correlation between the cluster evidence for pathway-pathway interactions and the Jaccard similarity that reflects number of genes in common between pathways (R-sqared ~ 0.02).
# 2. The pathway interactions with the highest cluster evidence mostly have small Jaccard similarity
# 3. Many small contributions add up to significant cluster pathway evidence. Discard filtering by individual clusters.
# 4. Filter by both Weight.clust and pathway Jaccard similarity
# ---------------------------------------------------------------------------------------------------

# Focus on cluster evidence but also graph bp edges
# Used normalized version on same scale (tpnn)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
filter.edges.between("EGF/EGFR signaling pathway","Transmembrane transport of small molecules", tpnn)
# Now 0.09122693  cluster evidence
filter.edges.between("EGF/EGFR signaling pathway","Glycolysis and gluconeogenesis", tpnn)
# Now 0.05601639 cluster evidence
filter.edges.0(c("EGF/EGFR signaling pathway","Glycolysis and gluconeogenesis", "Transmembrane transport of small molecules"), tpnn)
focus <- c("Glycolysis and gluconeogenesis", "EGF/EGFR signaling pathway", "Transmembrane transport of small molecules")
focus.df <- data.frame(id=focus)
focus.edges <- filter.edges.0(focus, pcnn)
intersect(bioplanet[["Glycolysis and gluconeogenesis"]], bioplanet[["Transmembrane transport of small molecules"]])     
#  11 genes pathway Jaccard similarity 0.02217742

focus.suid <- createNetworkFromDataFrames.check(focus.df, focus.edges, title="Focus Trio Normalized", collection = "Interactions")
setEdgeWidths.RCy32(focus.edges, factor=30, log=F)         
# BP edge very thin! with non-normalized, better with normalized, but linear is clearer:
setEdgeWidths.RCy32(focus.edges, factor=2, log=T)  
setEdgeSelectionColorDefault (col2hex("chartreuse"))
setDefaultNodeColor("#33FFFF" ) # bright Cyan
edgeColors <- c(col2hex(alpha("purple",0.33333)), col2hex("green"))
edgeColors <- c("#9966FF", col2hex("green"))
edgeTypes <- c("cluster evidence", "pathway Jaccard similarity")
setEdgeColorMapping( 'interaction', edgeTypes, edgeColors, 'd', default.color="#FFFFFF")
style.name <- "PCN style"
copyVisualStyle('default', style.name)
setVisualStyle(style.name)
# Expand using filtered networks
test1 <- filter.edges.1(focus, tpnn.no.bp)
# 2745
test2 <- filter.edges.1(focus, tpnn)
# 4041
# Filter for high cluster weight, low or high bp weight
pcnn.clust.hi <- pcnn.clust[pcnn.clust$Weight > 0.065,]
pcnn.bp.lo <- pcnn.bp[pcnn.bp$Weight < max(pcnn.bp$Weight)/10,] # lower than 0.07777778
pcnn.bp.hi <- pcnn.bp[pcnn.bp$Weight > 0.2,] # 3732
pcnn.bp.v.hi <- pcnn.bp[pcnn.bp$Weight > 0.5,] # 289
pcnn.filter1 <- rbind(pcnn.clust.hi, pcnn.bp.lo)
test3 <- filter.edges.1(focus, pcnn.clust.hi) # 141
test4 <- filter.edges.1(focus, pcnn.bp.hi) # 9
#focus.deg1.df <- data.frame(id=extract.gene.names.RCy3(test3)) # 73
focus.deg1.edges <- rbind(test3, test4)
focus.deg1.df <- data.frame(id=unique(c(focus.deg1.edges$source, focus.deg1.edges$target)))

focus.deg1.suid <- createNetworkFromDataFrames.check(focus.deg1.df, focus.deg1.edges, title="Focus Trio Degree 1", collection = "Interactions")
setEdgeWidths.RCy32(focus.edges, factor=2, log=T)  
# Filter edges
hist(focus.deg1.edges$Weight, col="red", breaks=40)
focus.deg1.bp <- focus.deg1.edges[which(focus.deg1.edges$interaction=="pathway Jaccard similarity"),]
# gets highly related pathways
focus.deg1.hiclust <- test3[which(test3$Weight>0.15),]
focus.deg.hiedges <- rbind(focus.edges, focus.deg1.bp, focus.deg1.hiclust)
focus.deg1.df2 <- data.frame(id=unique(c(focus.deg.hiedges$source, focus.deg.hiedges$target)))
focus.deg1.suid2 <- createNetworkFromDataFrames.check(focus.deg1.df2, focus.deg.hiedges, title="Focus Trio Degree 1, filtered", collection = "Interactions")


intersect(bioplanet[["EGFR1 pathway"]], bioplanet[["EGF/EGFR signaling pathway"]])
intersect(bioplanet[["Epidermal growth factor receptor (EGFR) pathway"]], bioplanet[["EGF/EGFR signaling pathway"]])
filter.edges.0(c("EGFR1 pathway","Epidermal growth factor receptor (EGFR) pathway","EGF/EGFR signaling pathway"), pcnn)
# cluster evidence range 0.03 to 0.073; bp 0.015 to 0.21
# How to represent the entire network for supplementary figure?
hist(pcnn$Weight, breaks=1000, col="gray88", ylim=c(0, 100))
hist(tpnn.no.bp$Weight, breaks=200, col="gray88", ylim=c(0, 100))
hist(pcnn.bp.hi$Weight, breaks=100, col="green")
tpnn.no.bp.hi <- tpnn.no.bp[tpnn.no.bp$Weight>0.05,] # 997
filtered.pcn <- rbind(tpnn.no.bp.hi, pcnn.bp.v.hi)
filtered.pcn.df <- data.frame(id=unique(c(filtered.pcn$source, filtered.pcn$target)))
focus.deg1.suid3 <- createNetworkFromDataFrames.check(filtered.pcn.df, filtered.pcn, title="PCN filtered", collection = "Interactions")
setEdgeWidths.RCy32(filtered.pcn, factor=2, log=T)  
# *** PCN_Networks_Focus.cys
# High BP edges are not unexpected and maybe distracting
tpnn.no.bp.hi.df <- data.frame(id=c(tpnn.no.bp.hi$source, tpnn.no.bp.hi$target))
focus.deg1.suid4 <- createNetworkFromDataFrames.check(tpnn.no.bp.hi.df, tpnn.no.bp.hi, title="PCN filtered, no bioplanet", collection = "Interactions")
setEdgeWidths.RCy32(tpnn.no.bp.hi, factor=30, log=F)         
setEdgeWidths.RCy32(tpnn.no.bp.hi, factor=2, log=T)  
# **** yFiles Organic layout with adjustments for figure

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---







tpn <- total.pathway.net
tpn.v.small <- tpn[tpn$Combined.Weight>10,] # 26
tpn.v.small[order(tpn.v.small$Weight.bp, decreasing=T),]
# mRNA processing, gene expression, translation, disease: big related clusters
tpn.small <- tpn[tpn$Combined.Weight>2,]
# 1389 edges <- 
filter.edges.1("EGF/EGFR signaling pathway", tpn.small)
tpn.small.nbp <-  tpn.nbp[tpn.nbp$Weight>2,]
# 84 edges

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
# Weight 1.66 (cluster evidence)
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
       
#****
# Now get cfn for each as above
look4 <- filter.edges.between( bioplanet[["Transmembrane transport of small molecules"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
# 13 edges
look5 <- filter.edges.between( bioplanet[["Glycolysis and gluconeogenesis"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
# 22 edges
# selectNodes(nodes=extract.gene.names.RCy3(look4), by.col="id")
# selectFirstNeighbors()
# createSubnetwork(subnetwork.name="Edges Between Pathays")
# Try this (again, in fresh Cytoscape session, with checking function)
focus.cccn1 <- graph.cfn.cccn.check (look4, ld=FALSE, gz=TRUE, only.cfn=FALSE)
focus.cccn2 <- graph.cfn.cccn.check (look5, ld=FALSE, gz=TRUE, only.cfn=FALSE)
# Merge networks in cytoscape, then
all.ratio.styles()
toggleGraphicsDetails()

#* PCN_WC2_WBPzeropoint3.cys, EGFRsignalingGlycolysisCFNCCCN_New.cys, EGFRsignalingTransportersCFNCCCN_New.cys
#___________________________________________________________________________
# What are the CFN/CCCN links between pathways linked by cluster evidence alone?
# gzallt.physical.network =	edge file of the physical interaction network: CCCN, CFN, plus negative corr edges
#_________________________________________
# Interactions with EGF/EGFR signaling pathway
# Choose Transmembrane transport of small molecules (druggable?)
a <- bioplanet[["Transmembrane transport of small molecules"]]
b <- bioplanet[["EGF/EGFR signaling pathway"]]
ab <- intersect (a, b) #0
ab.all <- unique(c(a, b)) # 573
# Method 1: composite shortest paths
gzpathgenes <- ab.all[ab.all %in% gzallt.gene.key$Gene.Name] # 201

egfr.transporters <- composite.shortest.paths(genes1=a[a %in% gzallt.gene.key$Gene.Name], genes2=b[b %in% gzallt.gene.key$Gene.Name], network=gzalltgene.physical.cfn.merged, exclude="MYH9")
cccn2 <- graph.cfn.cccn.check (egfr.transporters, ld=FALSE, gz=TRUE, only.cfn=FALSE)
# update to use drug medians
cccn3 <- graph.cfn.cccn.medians.check (egfr.transporters, ld=FALSE, gz=TRUE, only.cfn=FALSE)
# NOte: medians has many zeros (FIXED). Repeate using means
cccn4 <- graph.cfn.cccn.means.check (egfr.transporters, ld=FALSE, gz=TRUE, only.cfn=FALSE)
nodeDprops.RCy32()
FixEdgeDprops.RCy32()
#all.ratio.styles()
drug.medians.ratio.styles()
toggleGraphicsDetails()
# Save median styles file in EGF transporters networks 3.cys (medians) and 4 (means)
# look at CFN version
selectNodes(nodes=extract.gene.names(egfr.transporters), by.col="id", preserve.current.selection = F)
createSubnetwork(subnetwork.name="EGFR transporters composite shortest paths CFN")
# separate out
selectNodes(nodes=a, by.col="id", preserve.current.selection = F)
selectNodes(nodes=b, by.col="id", preserve.current.selection = F)
# What are the edges *between* uniuqe pathway genes?
look4 <- filter.edges.between( bioplanet[["Transmembrane transport of small molecules"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
# 13 edges
selectNodes(nodes=extract.gene.names.RCy3(look4), by.col="id")
selectFirstNeighbors()
createSubnetwork(subnetwork.name="Edges Between Pathways")
# Still to complex, try this (again)
cccn4 <- graph.cfn.cccn (look4, ld=FALSE, gz=TRUE, only.cfn=FALSE)
cccn5 <- graph.cfn.cccn.medians.check (look4, ld=FALSE, gz=TRUE, only.cfn=FALSE)
cccn6 <- graph.cfn.cccn.means.check (look4, ld=FALSE, gz=TRUE, only.cfn=FALSE)
setCorrEdgeAppearance(cccn4) 
setCorrEdgeAppearance(getTableColumns('edge'))
FixEdgeDprops.RCy32()
# Already done all.ratio.styles()
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
mutual.friends <- intersect(extract.gene.names(transporter.friends), extract.gene.names(EGFR.friends))
mutual.friends <- mutual.friends %w/o% c(a,b)
# 120 mutual friends 
# Network of these
fof.cfn <- filter.edges.0(c(a,b,friends.of.friends),gzalltgene.cfn.rcy3)
# 2634 edges
cfn1 <- graph.cfn.cccn (fof.cfn, ld=FALSE, gz=TRUE, only.cfn=TRUE)
toggleGraphicsDetails()
# >>>>>>>
# Use simpler merged version and filter.edges.between
transporter.friends.physical <- filter.edges.1(a, gzalltgene.physical.cfn.merged) 
EGFR.friends.physical <- filter.edges.1(b, gzalltgene.physical.cfn.merged) 
mutual.friends.physical <- intersect(extract.gene.names(transporter.friends.physical), extract.gene.names(EGFR.friends.physical))
mutual.friends.physical <- mutual.friends.physical %w/o% c(a,b)
# 54 genes
fof.cfn2 <- filter.edges.between(c(a,b),mutual.friends.physical, gzalltgene.physical.cfn.merged) # now 232 edges
cfn2 <- graph.cfn.cccn.check (fof.cfn2, ld=FALSE, gz=TRUE, only.cfn=TRUE)
#>>>>>
cfn3 <- graph.cfn.cccn.medians.check(fof.cfn2, ld=FALSE, gz=TRUE, only.cfn=TRUE)
cfn4 <- graph.cfn.cccn.means.check(fof.cfn2, ld=FALSE, gz=TRUE, only.cfn=TRUE)
nodeDprops.RCy32()
FixEdgeDprops.RCy32()
# separate out
selectNodes(nodes=a, by.col="id", preserve.current.selection = F)
selectNodes(nodes=b, by.col="id", preserve.current.selection = F)
# Method to determine whether mutual friends connect to both pathways
mutual.friends.cfn1 <- filter.edges.between(a, friends.of.friends, gzalltgene.physical.cfn.merged)
mutual.friends.cfn2 <- filter.edges.between(b, friends.of.friends, gzalltgene.physical.cfn.merged)
mutual.friends.cfn <- rbind(mutual.friends.cfn1, mutual.friends.cfn2)
# also 296 edges
identical(mutual.friends.cfn, fof.cfn2) # F
# Still has single friends
cfn3 <- graph.cfn.cccn.check (mutual.friends.cfn, ld=FALSE, gz=TRUE, only.cfn=TRUE)
# separate out
selectNodes(nodes=a, by.col="id", preserve.current.selection = F)
selectNodes(nodes=b, by.col="id", preserve.current.selection = F)
# For CFN/CCCN figure
setNodeLabelColorDefault(col2hex("transparent"), style.name = "all.drug Style",)

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
# *** EGF transporters networks.cys
# Who is where? And in what family
between.genes <- e.t.b.df[e.t.b.df$betweenness>0, "Gene.Name"]
intersect(between.genes,   bioplanet[["Transmembrane transport of small molecules"]]) # 12
intersect(between.genes,   bioplanet[["EGF/EGFR signaling pathway"]]) # 12
# 37; EGFR is at the top. Are these related to overall betweenness?
intersect(between.genes, mutual.friends)
#**** 77/120
intersect(between.genes, mutual.friends.physical) -> mutual.between
#****46/54
# 
selectNodes(nodes=mutual.between,  by.col="id", preserve.current.selection = F)
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
etb.df[unlist(sapply(etb.df, function(x) grep("RNA", x, fixed = T))),]
etb.RNA <- etb.df[unlist(sapply(etb.df, function(x) grep("RNA", x, fixed = T))), "Gene.Name"]
etb.df[unlist(sapply(etb.df, function(x) grep("kinase", x, fixed = T))),]
etb.kinase <- etb.df[unlist(sapply(etb.df, function(x) grep("kinase", x, fixed = T))), "Gene.Name"]
etb.df[unlist(sapply(etb.df, function(x) grep("biquitin", x, fixed = T))),]
etb.ubn <- etb.df[unlist(sapply(etb.df, function(x) grep("biquitin", x, fixed = T))), "Gene.Name"]
# Which are mutual friends?
mutual.friends.physical[mutual.friends.physical %in% etb.kinase]
mutual.friends.physical[mutual.friends.physical %in% etb.RNA]
mutual.friends.physical[mutual.friends.physical %in% etb.ubn]
etb.df[etb.df$Gene.Name %in% mutual.friends.physical[mutual.friends.physical %in% etb.ubn],]
etb.mfs <- etb.df[etb.df$Gene.Name %in% mutual.friends.physical,]
# NOTE: RAN is returned because grep finds anthing with the ame letters. 
selectNodes(mutual.friends.physical[mutual.friends.physical %in% etb.kinase])
  # For figures with no node labels
setNodeLabelColorDefault(col2hex("transparent"), style.name = "H1781_AfatinibRatio Style",)
#____________________________________________________________________________________________________________
#____________________________________________________________________________________________________________
# Interactions with EGF/EGFR signaling pathway
# Choose Glycolysis and gluconeogenesis (druggable?)
a <- bioplanet[["Glycolysis and gluconeogenesis"]]
b <- bioplanet[["EGF/EGFR signaling pathway"]]
ab <- intersect (a, b) #0
ab.all <- unique(c(a, b)) # 216
# Method 1: composite shortest paths
gzpathgenes <- ab.all[ab.all %in% gzallt.gene.key$Gene.Name] # 111
egfr.glycolosis <- composite.shortest.paths(genes1=a[a %in% gzallt.gene.key$Gene.Name], genes2=b[b %in% gzallt.gene.key$Gene.Name], network=gzalltgene.physical.cfn.merged, exclude="MYH9") #1431
# cccn5 <- graph.cfn.cccn.check (egfr.glycolosis, ld=FALSE, gz=TRUE, only.cfn=FALSE)
# Use medians
cccn5 <- graph.cfn.cccn.medians.check(egfr.glycolosis, ld=FALSE, gz=TRUE, only.cfn=FALSE)
nodeDprops.RCy32()
FixEdgeDprops.RCy32()
# all.ratio.styles()
drug.medians.ratio.styles()
toggleGraphicsDetails()
# look at CFN version
selectNodes(nodes=extract.gene.names(egfr.glycolosis), by.col="id", preserve.current.selection = F)
createSubnetwork(subnetwork.name="EGFR glycolosis composite shortest paths CFN")
# separate out
selectNodes(nodes=a, by.col="id", preserve.current.selection = F)
selectNodes(nodes=b, by.col="id", preserve.current.selection = F)

# What are the edges *between* uniuqe pathway genes?
look5 <- filter.edges.between( bioplanet[["Glycolysis and gluconeogenesis"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
# 22 edges
selectNodes(nodes=extract.gene.names.RCy3(look5), by.col="id")
selectFirstNeighbors()
createSubnetwork(subnetwork.name="Edges Between Pathways")
# Still to complex, try this (again)
cccn6 <- graph.cfn.cccn.check (look5, ld=FALSE, gz=TRUE, only.cfn=FALSE)
#setCorrEdgeAppearance(cccn6) 
#setCorrEdgeAppearance(getTableColumns('edge'))
FixEdgeDprops.RCy32()
# Already done all.ratio.styles()
# Who is in which pathway?
look5.df <- data.frame(Gene.Name=sort(unique(c(look5$source, look5$target))))
look5.df$EGFR.path <- sapply(look5.df$Gene.Name, function (x) x %in% b)
look5.df$Transporter.path <- sapply(look5.df$Gene.Name, function (x) x %in% a)
# >>>>>>>>>>>>>>>>>
# "Friends of friends" approach
# Use CFN here; keep genetic interactions (gzalltgene.cfn.rcy3)
a <- bioplanet[["Glycolysis and gluconeogenesis"]]
b <- bioplanet[["EGF/EGFR signaling pathway"]]
# glycolosis.friends <- filter.edges.1(a, gzalltgene.cfn.rcy3) 
# EGFR.friends <- filter.edges.1(b, gzalltgene.cfn.rcy3) 
# 142 mutual friends
# friends.of.friends <- intersect(extract.gene.names(glycolosis.friends), extract.gene.names(EGFR.friends)) %w/o% c(a,b)
glycolosis.friends <- filter.edges.1(a, gzalltgene.physical.cfn.merged) 
EGFR.friends <- filter.edges.1(b, gzalltgene.physical.cfn.merged) 
mutual.friends <- intersect(extract.gene.names(glycolosis.friends), extract.gene.names(EGFR.friends))
mutual.friends <- mutual.friends %w/o% c(a,b)
# 85 mutual friends
# Method to determine whether mutual friends connect to both pathways
mutual.friends.cfn1 <- filter.edges.between(a, mutual.friends, gzalltgene.physical.cfn.merged)
mutual.friends.cfn2 <- filter.edges.between(b, mutual.friends, gzalltgene.physical.cfn.merged)
mutual.friends.cfn <- rbind(mutual.friends.cfn1, mutual.friends.cfn2)
# 396 vs 461 edges w/genetic

cfn7 <- graph.cfn.cccn.check (mutual.friends.cfn, ld=FALSE, gz=TRUE, only.cfn=TRUE)
# Try with median
cfn8 <- graph.cfn.cccn.medians(mutual.friends.cfn, ld=FALSE, gz=TRUE, only.cfn=TRUE, pruned=FALSE)
# separate out
selectNodes(nodes=a, by.col="id", preserve.current.selection = F)
selectNodes(nodes=b, by.col="id", preserve.current.selection = F)
FixEdgeDprops.RCy32()

# Network of these
fof.cfn <- filter.edges.0(c(a,b,friends.of.friends),gzalltgene.cfn.rcy3)
# 3184 edges
cfn3 <- graph.cfn.cccn.check (fof.cfn, ld=FALSE, gz=TRUE, only.cfn=TRUE)
toggleGraphicsDetails()
# Use simpler merged version and filter.edges.between
fof.cfn2 <- filter.edges.between(c(a,b),friends.of.friends, gzalltgene.physical.cfn.merged) # now 461 edges
cfn4 <- graph.cfn.cccn.check (fof.cfn2, ld=FALSE, gz=TRUE, only.cfn=TRUE)
# separate out
selectNodes(nodes=a, by.col="id", preserve.current.selection = F)
selectNodes(nodes=b, by.col="id", preserve.current.selection = F)
# ****
# Try this with different views, including composite shortest paths and physical between in "EGF glycolosis Network.cys"
selectNodes(nodes=friends.of.friends, by.col="id", preserve.current.selection = F)
# Who is not in friends.of.friends in composite shortest paths: egfr.glycolosis
csp.genes <- outersect(extract.gene.names(egfr.glycolosis), c(a, b, friends.of.friends)) # 308 genes§
selectNodes(nodes=csp.genes, by.col="id", preserve.current.selection = F)
# Back to csp for CCCN version
selectNodes(nodes=extract.gene.names(egfr.glycolosis), by.col="id", preserve.current.selection = F)
invertNodeSelection()
createSubnetwork(subnetwork.name="EGFR glycolosis composite shortest paths CCCN")
# *** EGF glycolosis networks.cys

#
# Idea: betweenness of nodes between pathways
egfr.glycolosis.g <- graph_from_data_frame(egfr.glycolosis)
e.g.between <- betweenness(egfr.glycolosis.g)
e.g.b.df <- data.frame(Gene.Name=names(e.g.between), betweenness=e.g.between)
e.g.b.df <- e.g.b.df[order(e.g.b.df$betweenness, decreasing=TRUE),]
# Who is where? And in what family
between.genes <- e.g.b.df[e.g.b.df$betweenness>0, "Gene.Name"]
intersect(between.genes,   bioplanet[["Glycolysis and gluconeogenesis"]]) # 13
intersect(between.genes,   bioplanet[["EGF/EGFR signaling pathway"]]) # 37
# 37; EGFR is at the top. Are these related to overall betweenness?
egb.df <- left_join(e.g.b.df, gzgene.cfn.netatts, by="Gene.Name")
gz.cf.ppibetween <- gz.cf[which(gz.cf$Node.ID=="gene"), c("Gene.Name", "ppibetween", "norm.ppibetween", "HPRD.Function")]
egb.df <- left_join(egb.df, gz.cf.ppibetween[, c("Gene.Name", "HPRD.Function")], by="Gene.Name")
egb.df$enzyme.class <- sapply(egb.df$Gene.Name, function (x) names(zymes.list[grep(x, zymes.list)]))

# Relationships to pathway crosstalk network betweenness from CFN and PPI networks
plot(egb.df$betweenness~egb.df$ppibetween, pch=19, col="seagreen")
summary(lm(egb.df$betweenness~egb.df$ppibetween))
# R-squared:  0.2233 
plot(egb.df$betweenness~egb.df$norm.ppibetween, pch=19, col="forestgreen")
# norm.ppibetween is CFN ppibetween/allppibetween
summary(lm(egb.df$betweenness~egb.df$norm.ppibetween))
# R-squared: 0.02476
plot(egb.df$betweenness~egb.df$allppibetween, pch=19, col="forestgreen")
summary(lm(egb.df$betweenness~egb.df$allppibetween))
# R-squared: 0.006232

# Which genes are high betweenness relative to total CFN?
#
egb.df[which(egb.df$betweenness>egb.df$allppibetween), "Gene.Name"]
# [1] "CTTN"  "PGAM1"
egb.df[which(egb.df$betweenness>egb.df$ppibetween), "Gene.Name"] -> highbetweeners
#"CTTN"  "GJA1"  "PGAM1" "OCLN"  "BCAR1" "RALA"  "KRAS"  "GFM1" 
# Vs. above EGFR/transporters [1] "CTTN"   "EFNB1"  "OCLN"   "PAK2"   "GJA1"   "BCAR1"  "PRKACA" "PIK3R1" "MAPK14" "EPHB4" "RALA"   "PLAA"   "IL6ST"  "PSAT1"  "KIT"   

selectNodes(highbetweeners, by.col="id", preserve.current.selection = F)
selectFirstNeighbors()
# Ask if these are enzymes
egb.df[egb.df$Gene.Name %in% highbetweeners, c("Gene.Name", "HPRD.Function", "enzyme.class")]
# Group by enzyme or function
egb.RNA <- egb.df[grep("RNA", egb.df[,c("HPRD.Function","enzyme.class")], fixed = T), "Gene.Name"]
egb.kinase <- egb.df[grep("kinase", egb.df[,c("HPRD.Function","enzyme.class")]), "Gene.Name"]

# For figures with no node labels
setNodeLabelColorDefault(col2hex("transparent"), style.name = "H1781_AfatinibRatio Style",)


#------::::::::>>>>>>>>>>
# Another alternitive look: Find peptide edges between pathways.
get.pep.egdes.between.pathways <- function(pathway1, pathway2, only_between=FALSE, pepclusterlist=eu.sp.sed.gzallt, cccnedges=gzallt.cccnplus, pathway.evidence=cluster.pathway.evidence, ev.threshold=0){
  genes1 <- bioplanet[[pathway1]]
  genes2 <- bioplanet[[pathway2]]
  pathevidence.subset <- pathway.evidence[,c(pathway1, pathway2)]
  ev.list <- apply(pathevidence.subset, 1, function(y) y %>% keep(., function(x) !is.na(x) & x > ev.threshold) )
  if(length(ev.list) < 1) {stop("Threshold too high!")}
  # Look for more than two pathways
  ev.list.sizes <- sapply(ev.list, length)
  ev.list2 <- ev.list[which(ev.list.sizes>=2)]
  if(length(ev.list2) < 1) {stop("Threshold too high!")}
  # Clusters are 
  clusternames <- names(ev.list2)
  cluster.genes.list <- lapply(pepclusterlist[clusternames], get.gene.names.from.peps)
  a <- lapply(cluster.genes.list, function(x) intersect(x, genes1))
  b <- lapply(cluster.genes.list, function(x) intersect(x, genes2))
  # ab <- mapply(c, a, b)
  pepedges.list <- list()
  for (i in 1:length(clusternames)) {
    a.genes <- a[[i]] 
    b.genes <- b[[i]]
    clustpeps <- pepclusterlist[[clusternames[i]]]
    # Get all peptides from cluster
    a.peps <- NULL
    b.peps <- NULL
    for (j in 1:length(a.genes)) {
      a.peps[[j]] <- clustpeps[grep(a.genes[[j]], clustpeps)] }
    for (k in 1:length(b.genes)) {
      b.peps[[k]] <- clustpeps[grep(b.genes[[k]], clustpeps)] }
    if (only_between==TRUE) {
      pepedges.list[[i]] <- filter.edges.between(a.peps, b.peps, cccnedges)
    } else {
      pepedges.list[[i]] <- filter.edges.0(c(a.peps, b.peps), cccnedges)
    }
  }
  # Some edges are NA
  pepedges.list <- pepedges.list[!is.na(pepedges.list)]
  pepedges <- do.call("rbind", pepedges.list)
  return(pepedges)
}
tryit <- get.pep.egdes.between.pathways(pathway1="Transmembrane transport of small molecules", pathway2 ="EGF/EGFR signaling pathway")
# make gene-peptide edges
# onlynecessary if gzallt.cccn.edges.plus, instead of gzallt.cccnplus, is used:
# transegf.tryit <- edgeType.to.interaction(tryit)
genenames <- extract.gene.names(tryit) # alternatively:
genenames <- unique(sapply(c(tryit[,1],tryit[,2]),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
# NOTE: Find gz.cf unpruned should be 12K rows! Find on google drive LD_GZ Networks
# TenCell.RData - fixed and saved
cccn.cf <- gz.cf[gz.cf$Gene.Name %in% genenames,]
# Revised using drug medians, means
cccn.cf <- gz.drug.medians.pruned[gz.drug.medians.pruned$Gene.Name %in% genenames,]
cccn.cf <- gz.drug.means.pruned[gz.drug.means.pruned$Gene.Name %in% genenames,]

net.gpe <- data.frame(source=cccn.cf$Gene.Name, target=cccn.cf$id, Weight=0.25, interaction="peptide")
net.gpe <- remove.autophos.RCy3(net.gpe)
# Use only those edges above
transegf.tryit <- tryit
net.gpe.pruned <- net.gpe[net.gpe$target %in% unique(c(transegf.tryit$source, transegf.tryit$target)),]
transegf.between <-filter.edges.between( bioplanet[["Transmembrane transport of small molecules"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
transegf.fe0 <- filter.edges.0( c(bioplanet[["Transmembrane transport of small molecules"]], bioplanet[["EGF/EGFR signaling pathway"]]), edge.file=gzalltgene.physical.cfn.merged)
transegf.edges <- rbind(net.gpe.pruned, transegf.between, transegf.tryit)
transegf.edges2 <- rbind(net.gpe.pruned, transegf.fe0, transegf.tryit)
transegf.cf <- gz.cf[gz.cf$id  %in% unique(c(transegf.edges$source, transegf.edges$target)),]
transegf.cf2 <- gz.cf[gz.cf$id  %in% unique(c(transegf.edges2$source, transegf.edges2$target)),]
# Using medians, means:
transegf.cf <- gz.drug.medians[gz.drug.medians$id  %in% unique(c(transegf.edges$source, transegf.edges$target)),]
transegf.cf2 <- gz.drug.medians[gz.drug.medians$id  %in% unique(c(transegf.edges2$source, transegf.edges2$target)),]
transegf.cf <- gz.drug.means[gz.drug.means$id  %in% unique(c(transegf.edges$source, transegf.edges$target)),]
transegf.cf2 <- gz.drug.means[gz.drug.means$id  %in% unique(c(transegf.edges2$source, transegf.edges2$target)),]

# Check
outersect(transegf.cf$id, unique(c(transegf.edges$source, transegf.edges$target)))
outersect(transegf.cf2$id, unique(c(transegf.edges2$source, transegf.edges2$target)))
intersect(mutual.friends, extract.gene.names.RCy3(transegf.edges)) # Only "RAC2
# Note: not zero if using "pruned"
tryit.graph <- createNetworkFromDataFrames(transegf.cf, transegf.edges, "EGF transport pep edges between")
setNodeMapping(transegf.cf)
nodeDprops.RCy32()
edgeDprops.RCy32()
setCorrEdgeAppearance(transegf.edges) 
drug.medians.ratio.styles()
drug.means.ratio.styles()
toggleGraphicsDetails() 
setEdgeWidths.RCy32(getAllEdges(), factor=1.25)
# Note:there are few directly connected nodes from CFN made by searching for "gene" and making a subnetwork
selectAllEdges()
selectNodesConnectedBySelectedEdges()
connected <- getSelectedNodes()
# return to EGF transport pep edges between
selectNodes(nodes=connected, by="id", preserve.current.selection = F)
# Examine one cluster at a time
clusternodes <- getSelectedNodes()
clustergenes <- extract.genes.from.peplist(clusternodes)
# return to bigger network
selectNodes(nodes=c(clusternodes, clustergenes), by.col="id", preserve.current.selection = F)
 # createSubnetwork(nodes=c(clusternodes, clustergenes), subnetwork.name = "Cluster 1 CFN CCCN")
 # Doesn't work, do manually.

#-----
# Degree 1 next for comparison
tryit.graph2 <- createNetworkFromDataFrames.check(transegf.cf2, transegf.edges2, "EGF transport degree 1")
setdiff (transegf.cf$id, transegf.cf2$id)  # not identical
# first is a subset so don't need to
# Merge networks in cytoscape, then
nodeDprops.RCy32()
edgeDprops.RCy32()
setNodeMapping(transegf.cf2)
setCorrEdgeAppearance(transegf.edges2) 
FixEdgeDprops.RCy32()
toggleGraphicsDetails() 
setEdgeWidths.RCy32(getAllEdges(), factor=1.25)
# If not done: all.ratio.styles() or drug.means.ratio.styles()
# Conclusion: a different, possibly better view of the CFN/CCCN.
# Compare
#
focus.cccn1 <- graph.cfn.cccn.check (transegf.between, ld=FALSE, gz=TRUE, only.cfn=FALSE)
focus.nodes <- getAllNodes() # 126 vs. in between network above
between.nodes <- getAllNodes() # 243 nodes
outersect(focus.nodes, between.nodes)
# Getting edges between pathway peptides is superior! ****
# Doesn't return superfluous peptides; can plot CFN, CCCN, or both
# Save in "EGF transporters networks 2. cys"
# Arrange CFN/CCCN
a <- bioplanet[["Transmembrane transport of small molecules"]]
b <- bioplanet[["EGF/EGFR signaling pathway"]]
selectNodes(nodes=a, by.col="id", preserve.current.selection = F)
selectNodes(nodes=b, by.col="id", preserve.current.selection = F)
# What affecgts SLC transporters?
lookSLC <- gz.cf.pruned[grep("SLC", gz.cf.pruned$Gene.Name),]
#
# Test with glycolysis and gluconeogenesis...
# >>>>>>>>>>>>••••••••••••
tryit2 <- get.pep.egdes.between.pathways(pathway1="Glycolysis and gluconeogenesis", pathway2 ="EGF/EGFR signaling pathway") #101 edges
tryit2b <- get.pep.egdes.between.pathways(pathway1="Glycolysis and gluconeogenesis", pathway2 ="EGF/EGFR signaling pathway", only_between = TRUE) #49 edges
# Create a new function to make gene-peptide edges from a peptide edgefile
make.genepep.edges <- function(peptide.edgefile, cccn.cf=gz.cf) {
  peptides <- unique(c(peptide.edgefile$source, peptide.edgefile$target))
  genenames <- sapply(peptides,  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])
  net.gpe <- data.frame(source=genenames, target=peptides, Weight=0.25, interaction="peptide")
  return(net.gpe)
}
tryit2t.gpe <- make.genepep.edges(tryit2)
glucegf.between <-filter.edges.between( bioplanet[["Glycolysis and gluconeogenesis"]], bioplanet[["EGF/EGFR signaling pathway"]], edge.file=gzalltgene.physical.cfn.merged)
glucegf.fe0 <- filter.edges.0( c(bioplanet[["Glycolysis and gluconeogenesis"]], bioplanet[["EGF/EGFR signaling pathway"]]), edge.file=gzalltgene.physical.cfn.merged)
glucegf.edges <- rbind(tryit2t.gpe, glucegf.between, tryit2)
glucegf.edges2 <- rbind(tryit2t.gpe, glucegf.fe0, tryit2)
#
glucegf.cf <- gz.cf[gz.cf$id  %in% unique(c(glucegf.edges$source, glucegf.edges$target)),]
glucegf.cf2 <- gz.cf[gz.cf$id  %in% unique(c(glucegf.edges2$source, glucegf.edges2$target)),]
# Try with medians
glucegf.cf <- gz.drug.medians[gz.drug.medians$id  %in% unique(c(glucegf.edges$source, glucegf.edges$target)),]
glucegf.cf2 <- gz.drug.medians[gz.drug.medians$id  %in% unique(c(glucegf.edges2$source, glucegf.edges2$target)),]
# 
# Check
outersect(glucegf.cf$id, unique(c(glucegf.edges$source, glucegf.edges$target)))
# Now zero
tryit2.graph <- createNetworkFromDataFrames.check(glucegf.cf, glucegf.edges)
setNodeMapping(glucegf.cf)
nodeDprops.RCy32()
edgeDprops.RCy32()
setCorrEdgeAppearance(glucegf.edges) 
drug.medians.ratio.styles()
# drug.means.ratio.styles()
toggleGraphicsDetails() 
# FE 0
tryit2.graph2 <- createNetworkFromDataFrames.check(glucegf.cf2, glucegf.edges2)
setdiff (glucegf.cf$id, glucegf.cf2$id)  
# first is a subset so don't need to
# Merge networks in cytoscape, then
nodeDprops.RCy32()
edgeDprops.RCy32()
setNodeMapping(glucegf.cf2)
setCorrEdgeAppearance(glucegf.edges2) 
# all.ratio.styles()
# drug.medians.ratio.styles()
# ids are the same
toggleGraphicsDetails() 
#**** EGF Glycolysis Network.cys; EGFR glycolysis medians.cys
#*#_______________________________________________________________
#*#
#*#xperiment with graphing clusters
#*Focus on eu.sp.sed.gzallt.data[["130.141.111"]]
#*# Note essgzallt.data <- eu.sp.sed.gzallt.data
clusterdata <- eu.sp.sed.gzallt.data[["130.141.111"]]
# Note this contains raw data too
# Limited log2 ratios can be found here:
gzp.data.log2 <- gz.cf.pruned[,17:41]
rownames(gzp.data.log2) <- gz.cf.pruned$Gene.Name
# Note this contains genes too! Start with gzdata.allt
gzratios.allt.log2 <- gzdata.allt[, grep("atio", names(gzdata.allt))]
gzratios.allt.linear <- 2^gzratios.allt.log2
clusterdata.log2 <- clust.data.from.vec(eu.sp.sed.gzallt[["130.141.111"]], gzratios.allt.log2)
clusterdata.linear <- clust.data.from.vec(eu.sp.sed.gzallt[["130.141.111"]], gzratios.allt.linear)
dev.new()
x <- graph.cluster.3c(clusterdata.log2)
dev.new()
y <- graph.clust6d(clusterdata.log2)
z <- graph.clust6d(clusterdata.linear)
# Log2 definately better
a <- graph.clust5d(clusterdata.log2)
# This version scales the graph (not recommended:
b <- graph.clust6(clusterdata.log2)
# Deal with NA in the way done for tsne dissimilarity
# Also resize windows 
xd <- graph.clust6d(clusterdata.log2)
xd.l <- graph.clust6d.l(clusterdata.log2)
# Smaller font for tall graphs
xd.la <- graph.clust6d.la (clusterdata.log2)
# Still rows incomplete!
#
# This version makes a pdf file:
# Note error figure margins too large!
xd.l2 <- graph.clust6d.l2 (clusterdata.log2[1:10, 1:10], "/Users/_mark_/Dropbox/_Work/R_/_LINCS/_KarenGuolin/130.141.111test.pdf")
#  version of pdf file maker that doesn't sort
xd.l3 <- graph.clust6d.l3 (clusterdata.log2, "/Users/_mark_/Dropbox/_Work/R_/_LINCS/_KarenGuolin/130.141.111test.pdf")
# Revised 
# # see https://sebastianraschka.com/Articles/heatmaps_in_r.html 
# See this explanation here: https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html
 # https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
# https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html
#33333333 Not working:
graph.clust6d.la2.bad <- function(cluster.data.index, dendro="both") {
  #  options          dendrogram = c("both","row","column","none"),
  custer.df <- data.frame(cluster.data.index)
  if (any(grepl("Gene.Name", colnames(cluster.df)))) {
    if (!identical(rownames(cluster.df), cluster.df$Gene.Name)) {
      rownames(cluster.df) <- cluster.df$Gene.Name    }
    cluster.df	<- cluster.df[, -grep("Gene.Name", colnames(cluster.df))] 
  }
  if (ncol(cluster.df) < 2) {
    cat ("\n","This is a single sample cluster!", "\n") 
    return (cluster.df) } else {
      cluster.m <- data.matrix(cluster.df)
      # optional:
      #rbyheatcolors <- colorRampPalette(colors=c('#0000FF',  '#FFFF00'), bias=0.25, space="rgb", interpolate = "spline")
      rbyheatcolors <- colorRampPalette(colors=c('#3333FF', '#FFFF00'), bias=0.25, space="rgb", interpolate = "linear")
      # royal blue to yellow
      palette(c(rbyheatcolors(500)))
      # Use dendrogram to order and modified dist function that makes NA 2*max(is.na(x))
      dev.new()  # can resize manually
      heatmap.2(cluster.m, dendrogram=dendro, trace="none", symbreaks=TRUE, na.color="black", labRow=rownames(cluster.m), labCol=colnames(cluster.m), Rowv=TRUE, Colv=TRUE, distfun=dist2, hclustfun=hclust, scale="none", col=palette(), colsep=NULL, rowsep=NULL, sepwidth=c(0,0), revC=FALSE, keysize=0.5, density.info='histogram', cexCol = 0.3 + 1/log10(ncol(cluster.m)), denscol="green", main=names(cluster.data.index)[1], margins=c(14, 10)) 
      return (cluster.m) }
}

# test diverging vs sequential pallett
clusterdata.log2.pos <- abs(cluster.m[1:10, 1:10])
clusterdata.log2.neg <- -clusterdata.log2.pos
dev.new()
look <- graph.clust6d.la2(clusterdata.log2.pos[1:10, 1:10]) # ?
look <- graph.clust6d.l(clusterdata.log2.pos[1:10, 1:10])
dev.new()
look <- graph.clust6d.l(clusterdata.log2.neg[1:10, 1:10])

graph.cluster.3c(clusterdata.log2.pos)
# scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))


# From Cuneyt
graph.name <- paste(comp_path, "Dropbox/_Work/R_/_LINCS/_KarenGuolin/gzall.edges", sep="")
graph.name2 <- paste(comp_path, "Dropbox/_Work/R_/_LINCS/_KarenGuolin/gld.edges", sep="")
gz.g <- read_graph(graph.name, format="ncol")
ld.g <- read_graph(graph.name2, format="ncol")
gz.edges <- data.frame(as_edgelist(gz.g))
ld.edges <- data.frame(as_edgelist(ld.g))
names(gz.edges) <- c("source", "target")
gz.edges$interaction <- "igraph"
gz.edges$Weight <- 0.2
ld.edges$interaction <- "igraph"
ld.edges$Weight <- 0.2

gz.cf <- data.frame(id=names(V(gz.g)))
ld.cf <- data.frame(id=names(V(ld.g)))

tryit.graph <- createNetworkFromDataFrames(gz.cf, gz.edges, "Test Cy from igraph")

