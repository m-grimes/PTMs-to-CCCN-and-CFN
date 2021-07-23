#
# Experiment with UMAP
# see https://pair-code.github.io/understanding-umap/
#  Mark Grimes
#   May 11, 2020
#______________________________________________________________________________________________________
#
# Load large data files: correlation matrices
# iMac Pro
load(file="/Users/_mark_/Archive/Terra_Byte/R_Archive_2/_LINCS/GZTenCellMatrices2.RData")
# New Locations:  
# "/Users/_mark_/Dropbox/ZBackup/R_Backup/R_Archive_2/_LINCS/GZTenCellMatrices2.RData"
# " /Volumes/.timemachine/A52C7FA9-94A9-4792-9298-2C9BD45A5DA4/2021-06-02-122159.backup/2021-06-02-122159.backup/iMac_Pro_HD\ -\ Data/Users/_mark_/Archive/Terra_Byte/R_Archive_2/_LI NCS/GZTenCellMatrices2.RData"
# For previous work:
source (" /Users/_mark_/Dropbox/_Work/R_/MG_packages.R")
load(file=paste(comp_path, "Dropbox/_Work/R_/_LINCS/_Export_for_Paper/LC_TMT_Nets.RData", sep=""))
# 
data_path <- paste(comp_path, "Dropbox/_Work/R_/_LINCS/_KarenGuolin/", sep="")
code_path <- paste(comp_path, "Dropbox/_Work/R_/PTMs-to-CCCN-and-CFN/", sep="")
load(file=paste(data_path, "KGFunDataObjects.RData", sep=""))
load(file=paste(data_path, "TenCell.RData", sep=""))
load(file=paste(data_path, "GZ_PPI_Networks2.RData", sep=""))
load(file=paste(data_path, "LD_NewCFNCCCN_Networks.RData", sep=""))
load(file=paste(data_path,"drug_effects.RData", sep=""))
#------------------------------------------------------------------------------------------------------
# NOTE: # error: even after terminal:  pip install umap-learn[plot]
#------------------------------------------------------------------------------------------------------
# Try UMAP on distance and dissimilarity matrices
# Allt = gzdata combined with tencell trimmed. 
# gzallt.dist = as.matrix (dist (gzdata.allt), method = "euclidean")
# dissimilarity.gzallt <- 1 - abs(gzallt.cor)
# SED: combinde Euclid and Spearman w/o taking absolute value.	
# gzallt.sed <- (gzallt.dist.1 + diss.gzallt.noabs)/2
# t-SNE clusters are in the list: eu.sp.sed.gzallt
# Try first with SED, default settings
gzallt.umap.sed <- umap(gzallt.sed, config = umap.defaults, method = c("naive", "umap-learn"))
head(gzallt.umap.sed)
# Format is compartmentalized into a kind of list
sapply(gzallt.umap.sed, class)
# Plotting is not trivial, see https://umap-learn.readthedocs.io/en/latest/plotting.html
# and https://www.r-bloggers.com/running-umap-for-data-visualisation-in-r/
plot(gzallt.umap.sed$layout)
# ggplot(data.frame(gzallt.umap.sed$layout)) - needs more configuration

# see umap.defaults: n_components can be set to 3
umap.defaults$n_components <- 3
gzallt.umap.sed3D <- umap(gzallt.sed, config = umap.defaults, method = "naive")
plot3d(gzallt.umap.sed3D$layout)
# Not very impressive. In both Cases there are three or four outlier
umap.defaults$n_components <- 2
#
# Will it work with NAs?
 # gzallt.umap <- umap(gzdata.allt, config = umap.defaults, method = c("naive", "umap-learn"))
   # No   " missing value where TRUE/FALSE needed"
# Try euclidean, since UMAP uses euclid
gzallt.umap.eu <- umap(gzallt.dist, config = umap.defaults, method = c("naive", "umap-learn"))
plot(gzallt.umap.eu$layout)
eu.umap <- data.frame(gzallt.umap.eu$layout)
# Which outlier?
gzallt.umap.eu$layout[which(gzallt.umap.eu$layout<(-100)),]
# Spearman dissimilarity
gzallt.umap.sp <- umap(dissimilarity.gzallt, config = umap.defaults, method = c("naive", "umap-learn"))
plot(gzallt.umap.sp$layout)
sp.umap <- data.frame(gzallt.umap.sp$layout)
sp.umap[which(sp.umap[,2]>60),]
# Group of 6
# EDNRA ubi K299, LMO7 p Y1667, RPL7 ubi K77, SRSF9 p Y179, SV2A ubi K143, TUBA1B ubi K394; TUBA4A ubi K394; TUBA1C ubi K394
gzdata.allt[rownames(sp.umap[which(sp.umap[,2]>60),]),]
# These PTMs are only in H3122 cells

# Rainbow plot of t-SNE clusters (as for ld data in LINCS_6.R)
plot(sp.umap)
for(i in 1:length(eu.sp.sed.gzallt))	{
  tp <- which(rownames(sp.umap) %in% eu.sp.sed.gzallt[[i]])
  points(sp.umap[tp,1], sp.umap[tp,2], col=rainbow(length(eu.sp.sed.gzallt))[i], pch=19)					}
# with 839 clusters the colors seem spread over many clusters
# Compare with t-SNE
head(sed.gzallt.tsne) # no rownames
# Zoom in for a closder look at dense areas
dev.new()
plot(sp.umap, ylim=c(-8, 18), xlim=c(-10, 10), col="lightgrey")
for(i in 1:length(eu.sp.sed.gzallt))	{
  if (i %in% seq(1, length(eu.sp.sed.gzallt), by=10)) plot(sp.umap, ylim=c(-8, 18), xlim=c(-10, 10), col="lightgrey")
  tp <- which(rownames(sp.umap) %in% eu.sp.sed.gzallt[[i]])
  points(sp.umap[tp,1], sp.umap[tp,2], col=sample(rainbow(length(eu.sp.sed.gzallt)),1), pch=19)					
  Sys.sleep(0.05) 
  if (i %in% seq(1, length(eu.sp.sed.gzallt), by=10)) Sys.sleep(1)  }
plot(eu.umap, ylim=c(-20, 20), col="lightgrey")
for(i in 1:length(eu.sp.sed.gzallt))	{
  tp <- which(rownames(eu.umap) %in% eu.sp.sed.gzallt[[i]])
  points(eu.umap[tp,1], eu.umap[tp,2], col=sample(rainbow(length(eu.sp.sed.gzallt)),1), pch=19)
  Sys.sleep(0.01)}
# This is interesting. Most t-SNE clusters are neighbors in UMAP plot, which is good. 
# Note that the minimum spanning tree, nearest neighbor clustering method would NOT work on these 2D plots.
# Is 3D useful?
sed3d.umap <- data.frame(gzallt.umap.sed3D$layout)
plot3d(sed3d.umap, ylim=c(-20, 20), col="lightgrey")

for(i in 1:length(eu.sp.sed.gzallt))	{
  tp <- which(rownames(eu.umap) %in% eu.sp.sed.gzallt[[i]])
  points3d(eu.umap[tp,1], eu.umap[tp,2], col=sample(rainbow(length(eu.sp.sed.gzallt)),1), pch=19)
  Sys.sleep(0.01)}


######################
# Are clusters which return similar pathays located in the same region on UMAP embeddings?
#

######################
# plot groups in 3D
load(file="/Users/_mark_/Archive/Terra_Byte/R_Archive_2/_LINCS/GZTenCellMatrices2.RData")
rownames(sed.gzallt.tsne) <- rownames(gzallt.dist)
sedcolors <- sample(rainbow(length(eu.sp.sed.gzallt)))
plot3d(sed.gzallt.tsne[rownames(sed.gzallt.tsne) %in% unlist(eu.sp.sed.gzallt),], col="grey88", axes=F, box=F, xlab="", ylab="", zlab="")
# Loop to color clusters
for (i in 1:length(eu.sp.sed.gzallt)) {
  #points3d(sed.tsne[eu.sp.sed.lfSYratios[[i]]], radius=0.8, col=sedcolors[i], add=T)
  tp <- which(rownames(sed.gzallt.tsne) %in% eu.sp.sed.gzallt[[i]])
  plot3d(sed.gzallt.tsne[tp,], type="s", radius=1, col=sedcolors[i], add=TRUE)		
}
# Took a screenshot.
# Note high resolution image can be saved
# https://stackoverflow.com/questions/38907514/saving-a-high-resolution-image-in-r