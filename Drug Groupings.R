#Hi Mark,Here are the groupings that I have been looking at. Also, I have been "averaging" using median rather than mean because I thought it would make the results less sensitive to a few extreme values. I added the all.drug group (#7) after talking to Eric Haura because he seemed interested in general responses to TKIs.
#  
#  Best,Karen
drug.groups <- list(  
  h3122.criz = c(
  "H3122SEPTM.C1.ratio",
  "H3122SEPTM.C2.ratio",
  "H3122SEPTM.C3.ratio",
  "H3122CrizotinibRatio",
  "H3122_Crizotinib.1Ratio"),
  hcc78.criz = c(
  "HCC78_CrizotinibRatio",
  "HCC78_Crizotinib.1Ratio"),
  pc9.erl = c(
  "PC9SEPTM.E1.ratio",
  "PC9SEPTM.E2.ratio",
  "PC9SEPTM.E3.ratio",
  "PC9_Erlotinib.1Ratio",
  "PC9_ErlotinibRatio"),
  all.criz = c(
  "H3122SEPTM.C1.ratio",
  "H3122SEPTM.C2.ratio",
  "H3122SEPTM.C3.ratio",
  "H3122CrizotinibRatio",
  "H3122_Crizotinib.1Ratio",
  "HCC78_CrizotinibRatio",
  "HCC78_Crizotinib.1Ratio",
  "H2228CrizotinibRatio",
  "STE.1_CrizotinibRatio"),
  all.erl = c(
  "PC9SEPTM.E1.ratio",
  "PC9SEPTM.E2.ratio",
  "PC9SEPTM.E3.ratio",
  "PC9_Erlotinib.1Ratio",
  "PC9_ErlotinibRatio",
  "HCC827_ErlotinibRatio",
  "HCC4006_ErlotinibRatio"),
  all.dasat = c(
  "H2286_DasatinibRatio",
  "H366_DasatinibRatio"),
  all.drug = c(
  "H3122SEPTM.C1.ratio",
  "H3122SEPTM.C2.ratio",
  "H3122SEPTM.C3.ratio",
  "H3122CrizotinibRatio",
  "H3122_Crizotinib.1Ratio",
  "HCC78_CrizotinibRatio",
  "HCC78_Crizotinib.1Ratio",
  "H2228CrizotinibRatio",
  "STE.1_CrizotinibRatio",
  "PC9SEPTM.E1.ratio",
  "PC9SEPTM.E2.ratio",
  "PC9SEPTM.E3.ratio",
  "PC9_Erlotinib.1Ratio",
  "PC9_ErlotinibRatio",
  "HCC827_ErlotinibRatio",
  "HCC4006_ErlotinibRatio",
  "H2286_DasatinibRatio",
  "H366_DasatinibRatio",
  "H1781_AfatinibRatio"))
 
# Use gz.cf to calculate medians
# Note: The values that are exactly 0 are really NA
# Convert to numeric class
options(stringsAsFactors=FALSE)
identical(names(gz.cf), names(gz.cf.pruned)) # TRUE
gz.cf.data <- data.matrix(gz.cf[, c(10:41, 44:46)])
gz.cf.head <- gz.cf[, -c(10:41, 44:46)]
gz.cf.NA <-  zero.to.NA(gz.cf.data)
gz.data.medians <- gz.cf.NA[, outersect(names(gz.cf.NA), unlist(drug.groups))]
for(i in 1:length(drug.groups)){
  gz.data.medians[,names(drug.groups)[i]] <- apply(gz.cf.NA[,drug.groups[[i]]], 1, function(x) median(x, na.rm=TRUE))
}
# Convert back into zeros for Cytoscape
gz.data.medians[is.na(gz.data.medians)] <- 0
# replace head
gz.drug.medians <- cbind(gz.cf.head, gz.drug.medians)
# Since this is for graphing also create pruned version that has only cccn nodes
gz.drug.medians.pruned <- gz.drug.medians[gz.drug.medians$id %in% gz.cf.pruned$id,]

# Funtion to make Cytoscape styles from this
drug.medians.ratio.styles <- function(ratiocols=NULL) {
  nodevalues <- getTableColumns('node')
  if(length(ratiocols)==0) {
    drugcolnames <- names(drug.medians)}
  for (i in 1:length(drug.medians)){ 
    plotcol <- drugcolnames[i]
    style.name = paste(drugcolnames[i], "Style")
    print(style.name)
    setVisualStyle("default")
    setNodeColorToRatios(plotcol)    
    copyVisualStyle('default', style.name)
    setVisualStyle(style.name)
  }
}
# Modify graphing function
graph.cfn.cccn.medians <- function(edgefile, ld=FALSE, gz=TRUE, only.cfn=FALSE, pruned=TRUE) {
  genenames <- extract.gene.names(edgefile)
  if(pruned==TRUE) {gz.cf <- gz.drug.medians.pruned} else {gz.cf <- gz.drug.medians}
  if (gz==TRUE) {
    cccn <- gzallt.cccnplus
    cccn.cf <- gz.cf[gz.cf$Gene.Name %in% genenames,]
  }
  if (ld==TRUE) {
    cccn <- ld.cccnplus
    cccn.cf <- ld.cf[ld.cf$Gene.Name %in% genenames,]
  }
  if (only.cfn==TRUE) {
    cfn.cf <- cccn.cf[which(cccn.cf$Node.ID=="gene"),]
    gene.suid <- createNetworkFromDataFrames(cfn.cf, edgefile, title=paste("CFN", (getNetworkCount()+1)), collection = "Interactions")
    setNodeMapping(cfn.cf)
    setCorrEdgeAppearance(edgefile)     
  }
  if (only.cfn==FALSE) {
    netpeps <- cccn.cf[which(cccn.cf$Node.ID=="peptide"), 'id']
    # make gene-peptide edges
    net.gpe <- data.frame(source=cccn.cf$Gene.Name, target=cccn.cf$id, Weight=0.25, interaction="peptide")
    # remove gene-gene interactions
    net.gpe <- remove.autophos.RCy3(net.gpe)
    ptm.cccn <-	filter.edges.0.RCy3(netpeps, cccn) 
    cfn.cccn.edges <- rbind(net.gpe, ptm.cccn, edgefile)
    if (gz==TRUE) {all.cf <- gz.cf[gz.cf$id  %in% unique(c(cfn.cccn.edges$source, cfn.cccn.edges$target)),]}
    if (ld==TRUE) {all.cf <- ld.cf[ld.cf$id  %in% unique(c(cfn.cccn.edges$source, cfn.cccn.edges$target)),]}
    cfn.cccn.suid <- createNetworkFromDataFrames(all.cf, cfn.cccn.edges, title=paste("CFN plus CCCN", (getNetworkCount()+1)), collection = "Interactions") 
    setNodeMapping(cccn.cf)
    nodeDprops.RCy32()
    setCorrEdgeAppearance(cfn.cccn.edges) 
    FixEdgeDprops.RCy32()
  }
  layoutNetwork("genemania-force-directed") 
  if (only.cfn==FALSE) return(cfn.cccn.edges)
}

graph.cfn.cccn.medians.check <- function(edgefile, ld=FALSE, gz=TRUE, only.cfn=FALSE, pruned=TRUE) {
  for (i in 1:10){
    graph.cfn.cccn.medians(edgefile, ld, gz, only.cfn, pruned) 
    edgeTable <- getTableColumns("edge", c("name", "shared name"))
    if(!identical(edgeTable[,1], edgeTable[,2])) {
      print(paste('Network', i, "is bad."))
      deleteNetwork()} else {
        print (paste("Network", i, "passes edge test."))
        break }
  }}
#
# Repeat to look at means
# Use gz.cf.NA to calculate means
gz.data.means <- gz.cf.NA[, outersect(names(gz.cf.NA), unlist(drug.groups))]
#
for(i in 1:length(drug.groups)){
  gz.data.means[,names(drug.groups)[i]] <- apply(gz.cf.NA[,drug.groups[[i]]], 1, function(x) mean(x, na.rm=TRUE))
}
# Since this is for graphing use pruned version that has only cccn nodes
# Convert back into zeros for Cytoscape
gz.data.means[is.na(gz.data.means)] <- 0
# replace head
gz.drug.means <- cbind(gz.cf.head, gz.data.means)
gz.drug.means.pruned <- gz.drug.means[gz.drug.means$id %in% gz.cf.pruned$id,]


# Compare means and medians
plot(unlist(gz.data.medians[,17:23]),  unlist(gz.data.means[,17:23]), pch=19, col=alpha("blue", 0.24))
# Funtion to make Cytoscape styles from this
drug.means.ratio.styles <- function(ratiocols=NULL) {
  nodevalues <- getTableColumns('node')
  if(length(ratiocols)==0) {
    drugcolnames <- names(drug.means)}
  for (i in 1:length(drug.means)){ 
    plotcol <- drugcolnames[i]
    style.name = paste(drugcolnames[i], "Style")
    print(style.name)
    setVisualStyle("default")
    setNodeColorToRatios(plotcol)    
    copyVisualStyle('default', style.name)
    setVisualStyle(style.name)
  }
}
# Modify graphing function
graph.cfn.cccn.means <- function(edgefile, ld=FALSE, gz=TRUE, only.cfn=FALSE, pruned=TRUE) {
  genenames <- extract.gene.names(edgefile)
  if(pruned==TRUE) {cf <- gz.drug.means.pruned} else {cf <- gz.drug.means}
  if (gz==TRUE) {
    cccn <- gzallt.cccnplus
    cccn.cf <- cf[cf$Gene.Name %in% genenames,]
  }
  if (ld==TRUE) { 
    cccn <- ld.cccnplus
    cccn.cf <- ld.cf[ld.cf$Gene.Name %in% genenames,]
  }
  if (only.cfn==TRUE) {
    cfn.cf <- cccn.cf[which(cccn.cf$Node.ID=="gene"),]
    gene.suid <- createNetworkFromDataFrames(cfn.cf, edgefile, title=paste("CFN", (getNetworkCount()+1)), collection = "Interactions")
    setNodeMapping(cfn.cf)
    setCorrEdgeAppearance(edgefile)     
  }
  if (only.cfn==FALSE) {
    netpeps <- cccn.cf[which(cccn.cf$Node.ID=="peptide"), 'id']
    # make gene-peptide edges
    net.gpe <- data.frame(source=cccn.cf$Gene.Name, target=cccn.cf$id, Weight=0.25, interaction="peptide")
    # remove gene-gene interactions
    net.gpe <- remove.autophos.RCy3(net.gpe)
    ptm.cccn <-	filter.edges.0.RCy3(netpeps, cccn) 
    cfn.cccn.edges <- rbind(net.gpe, ptm.cccn, edgefile)
    if (gz==TRUE) {all.cf <- cf[cf$id  %in% unique(c(cfn.cccn.edges$source, cfn.cccn.edges$target)),]}
    if (ld==TRUE) {all.cf <- ld.cf[ld.cf$id  %in% unique(c(cfn.cccn.edges$source, cfn.cccn.edges$target)),]}
    cfn.cccn.suid <- createNetworkFromDataFrames(all.cf, cfn.cccn.edges, title=paste("CFN plus CCCN", (getNetworkCount()+1)), collection = "Interactions") 
    setNodeMapping(cccn.cf)
    #nodeDprops.RCy32()
    setCorrEdgeAppearance(cfn.cccn.edges) 
    #FixEdgeDprops.RCy32()
  }
  layoutNetwork("genemania-force-directed") 
  if (only.cfn==FALSE) return(cfn.cccn.edges)
}

graph.cfn.cccn.means.check <- function(edgefile, ld=FALSE, gz=TRUE, only.cfn=FALSE, pruned=TRUE) {
  for (i in 1:10){
    graph.cfn.cccn.means(edgefile, ld, gz, only.cfn, pruned) 
    edgeTable <- getTableColumns("edge", c("name", "shared name"))
    if(!identical(edgeTable[,1], edgeTable[,2])) {
      print(paste('Network', i, "is bad."))
      deleteNetwork()} else {
        print (paste("Network", i, "passes edge test."))
        break }
  }}
