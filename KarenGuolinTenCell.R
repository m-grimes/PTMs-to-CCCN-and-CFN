#########################################################################################
# Pick up from KarenGuolinData4.R: read in ten cell line data
# 
load(file="/Users/mark_grimes/Software\ Images/GZTenCellMatrices3.RData") # Laptop
load(file="/Volumes/Terra_Byte/R_Archive_2/_LINCS/GZTenCellMatrices3.RData") # Office
load(file=paste("_LINCS/_KarenGuolin/", "TenCell.RData", sep=""))
load(file=paste("_LINCS/_KarenGuolin/", "GZ_PPI_Networks2.RData", sep=""))
#########################################################################################
  #   
tencellpath <- paste(comp_path, "Dropbox/_Work/R_/_LINCS/_KarenGuolin/10celllines_short_mapped_shorten_with_Localization0.8andPEP0.05cutoff_05292019",  sep="")
tencellphosname <- "TenCellPhosphoSites_mapped_shorten.txt"
tencellphos <- read.table(paste(tencellpath, tencellphosname, sep="/"), sep = "\t", skip = 0, header=TRUE, blank.lines.skip=T, fill=T, quote="\"", dec=".", comment.char = "", stringsAsFactors=F)
dim(tencellphos)
# 5941  383; after curation 3404
any(is.na(tencellphos$LeadingGeneSymbols))  # [1] TRUE
# there are some NAs in LeadingGeneSymbol
tail(tencellphos[order(tencellphos$LeadingGeneSymbols),  c(1:3, ncol(tencellphos))], 35 )
# Delete the NA rows
tencellphos <- tencellphos[!is.na(tencellphos$LeadingGeneSymbols),]
dim(tencellphos)  # 3373
headercols <- c("LeadingGeneSymbols", "Amino.acid", "Positions")
tencellphos.head <- tencellphos[,headercols]
# 
# Make peptide names using this function (revised to Karen's reformatting): 
name.peptide <- function (genes, modification="p", sites, aa)	{
  genes.v <- unlist(strsplit(genes, ";", fixed = TRUE))
  sites.v <- unlist(strsplit(sites, ";", fixed = TRUE))
  sites.v <- sapply(sites.v, function (x) paste (aa, x, sep=""))
  Peptide.v <- as.character(noquote(paste(genes.v[1:length(genes.v)], modification, sites.v[1:length(sites.v)], sep=" ")))
  Peptide <- paste(unique(Peptide.v), collapse="; ")
  return(Peptide)
}
tencellphos.head[is.na(tencellphos.head$LeadingGeneSymbols),]
# All gone
tencellphos.head[grep("-Sep", tencellphos.head$LeadingGeneSymbols),] # Some
tencellphos.head[grep("-Oct", tencellphos.head$LeadingGeneSymbols),] # 0
tencellphos.head[grep("-Mar", tencellphos.head$LeadingGeneSymbols),] # 0
# Excel goo: used fix.excel()
tencellphos.head$LeadingGeneSymbols <- sapply(tencellphos.head$LeadingGeneSymbols, fix.excel)
# Okay!
tencellphos.head$Peptide.Name <- mapply(name.peptide, genes=tencellphos.head$LeadingGeneSymbols, sites= tencellphos.head$Positions, aa=tencellphos.head$Amino.acid)
length(intersect(tencellphos.head$Peptide.Name, rownames(ld.ratio))) 
# 1144 ; mpw 1029 
length(intersect(tencellphos.head$Peptide.Name, rownames(kgp)))
# 676 ; now 589
# Check duplicates
any(duplicated(tencellphos.head$Peptide.Name)) # T
tencellphos.head$Peptide.Name[duplicated(tencellphos.head$Peptide.Name)]
# "HLA-A p Y344" "HLA-A p Y344"
tencellphos.head[grep("HLA-A", tencellphos.head$Peptide.Name), ]
tencellphos[grep("HLA-A", tencellphos$LeadingGeneSymbols), grep("Intensity", names(tencellphos))]
tencellphos.head[grep("HLA-A p Y344", tencellphos.head$Peptide.Name), ]
# These have different alpha chain isoforms that share a common peptide seq. Name by the number of possibilities
tencellphos.head[grep("HLA-A p Y344", tencellphos.head$Peptide.Name)[1], "Peptide.Name"] <- "HLA-A.4 p Y344"
tencellphos.head[grep("HLA-A p Y344", tencellphos.head$Peptide.Name)[2], "Peptide.Name"] <- "HLA-A.2 p Y344"
any(is.na(tencellphos.head$Peptide.Name))  # [1] FALSE
tencellphosdata <- tencellphos[,grep("Intensity", names(tencellphos))]
names(tencellphosdata) <- sapply(names(tencellphosdata), function (x) unlist(strsplit(x, "Intensity."))[2])
# To make column names identical, remove dates on the left and the last three characters
names(tencellphosdata) <- sapply(names(tencellphosdata), function (x) unlist(strsplit(x, "GZ_"))[2])
tcdnames <- names(tencellphosdata)
tcdnames <- make.unique(str_sub(tcdnames, 1, str_length(tcdnames)-3))
names(tencellphosdata) <- tcdnames
# make zero into NA, which it is.
zer0.tc <- which(tencellphosdata==0, arr.ind = TRUE)
tencellphosdata <- replace (tencellphosdata, zer0.tc, NA)
nmissing(tencellphosdata)/(dim(tencellphosdata)[1]*dim(tencellphosdata)[2])
# 71% NA
# Add rownames, first check above
# add rownames to data
rownames(tencellphosdata) <- tencellphos.head$Peptide.Name
hist(tencellphosdata, col="yellow", breaks=100)
########################################################################################
# Repeat with aceltylation data
tencellackname <- "AcetylSites_mapped_shorten_.txt"
tencellack <- read.table(paste(tencellpath, tencellackname, sep="/"), sep = "\t", skip = 0, header=TRUE, blank.lines.skip=T, fill=T, quote="\"", dec=".", comment.char = "", stringsAsFactors=F)
dim(tencellack)
# 2568  after curation 
# Delete the NA rows if necessary
tencellack[is.na(tencellack$LeadingGeneSymbols),] # 0
any(is.na(tencellack$LeadingGeneSymbols)) # F
headercols <- c("LeadingGeneSymbols", "Amino.acid", "Positions")
tencellack.head <- tencellack[,headercols]
# 
tencellack.head[is.na(tencellack.head$LeadingGeneSymbols),]
# All gone
tencellack.head[grep("-Sep", tencellack.head$LeadingGeneSymbols),] # Some
tencellack.head[grep("-Oct", tencellack.head$LeadingGeneSymbols),] # 0
tencellack.head[grep("-Mar", tencellack.head$LeadingGeneSymbols),] # 0
# Excel goo: used fix.excel()
tencellack.head$LeadingGeneSymbols <- sapply(tencellack.head$LeadingGeneSymbols, fix.excel)
# Okay!
tencellack.head$Peptide.Name <- mapply(name.peptide, genes=tencellack.head$LeadingGeneSymbols, sites= tencellack.head$Positions, modification="ack", aa=tencellack.head$Amino.acid)
length(intersect(tencellack.head$Peptide.Name, rownames(ld.ratio))) 
# 528 after fix.excel 
length(intersect(tencellack.head$Peptide.Name, rownames(kgdata)))
# 464
# Check duplicates
any(duplicated(tencellack.head$Peptide.Name)) # F
tencellack.head$Peptide.Name[duplicated(tencellack.head$Peptide.Name)] # 0
tencellackdata <- tencellack[,grep("Intensity", names(tencellack))]
dim(tencellackdata) # 2568   26
names(tencellackdata) <- sapply(names(tencellackdata), function (x) unlist(strsplit(x, "Intensity."))[2])

# make zero into NA, which it is.
zer0.tc <- which(tencellackdata==0, arr.ind = TRUE)
tencellackdata <- replace (tencellackdata, zer0.tc, NA)
nmissing(tencellackdata)/(dim(tencellackdata)[1]*dim(tencellackdata)[2])
# 73% NA
# To make column names identical, remove dates on the left and the last three characters
names(tencellackdata) <- sapply(names(tencellackdata), function (x) unlist(strsplit(x, "GZ_"))[2])
tcdnames <- names(tencellackdata)
tcdnames <- make.unique(str_sub(tcdnames, 1, str_length(tcdnames)-3))
names(tencellackdata) <- tcdnames
identical(names(tencellphosdata), names(tencellackdata))  # FALSE
data.frame(names(tencellphosdata), names(tencellackdata))
# Column names are out of order!
# add rownames to data
rownames(tencellackdata) <- tencellack.head$Peptide.Name
########################################################################################
# Repeat with ubiquitination data
tencellubname <- "GlyGlySites_mapped_shorten.txt"
tencellub <- read.table(paste(tencellpath, tencellubname, sep="/"), sep = "\t", skip = 0, header=TRUE, blank.lines.skip=T, fill=T, quote="\"", dec=".", comment.char = "", stringsAsFactors=F)
dim(tencellub) # 4736  383
# Delete the NA rows if necessary
any(is.na(tencellub$LeadingGeneSymbols))  # TRUE
tencellub <- tencellub[!is.na(tencellub$LeadingGeneSymbols),]
dim(tencellub)  # 4683
headercols <- c("LeadingGeneSymbols", "Amino.acid", "Positions")
tencellub.head <- tencellub[,headercols]
tencellub.head[is.na(tencellub.head$LeadingGeneSymbols),]
# All gone
tencellub.head[grep("-Sep", tencellub.head$LeadingGeneSymbols),] # Some
tencellub.head[grep("-Oct", tencellub.head$LeadingGeneSymbols),] # 0
tencellub.head[grep("-Mar", tencellub.head$LeadingGeneSymbols),] # 0
# Excel goo: used fix.excel()
tencellub.head$LeadingGeneSymbols <- sapply(tencellub.head$LeadingGeneSymbols, fix.excel)
# Okay!# 
tencellub.head$Peptide.Name <- mapply(name.peptide, genes=tencellub.head$LeadingGeneSymbols, sites= tencellub.head$Positions, modification="ubi", aa=tencellub.head$Amino.acid)
length(intersect(tencellub.head$Peptide.Name, rownames(ld.ratio))) 
# 0 
length(intersect(tencellub.head$Peptide.Name, rownames(kgdata)))
# 625
# Check duplicates
any(duplicated(tencellub.head$Peptide.Name)) # T
tencellub.head$Peptide.Name[duplicated(tencellub.head$Peptide.Name)]
# several HLAs:
# "HLA-A ubi K364" "HLA-A ubi K340" "HLA-A ubi K364" "HLA-C ubi K365" "HLA-A ubi K340" "HLA-B ubi K340"
#
tencellub.head[grep("HLA-A", tencellub.head$Peptide.Name), ]
# These have different alpha chain isoforms that share a common peptide seq. Name by the number of possibilities
tencellub.head[grep("HLA-A ubi K364", tencellub.head$Peptide.Name)[1], "Peptide.Name"] <- "HLA-A.4 ubi K364"
tencellub.head[grep("HLA-A ubi K364", tencellub.head$Peptide.Name)[1], "Peptide.Name"] <- "HLA-A.2 ubi K364"
tencellub.head[grep("HLA-A ubi K340", tencellub.head$Peptide.Name)[1], "Peptide.Name"] <- "HLA-A.4 ubi Y344"
tencellub.head[grep("HLA-A ubi K340", tencellub.head$Peptide.Name)[1], "Peptide.Name"] <- "HLA-A.2 ubi Y344"
any(duplicated(tencellub.head$Peptide.Name)) # T
tencellub.head[grep("HLA-B", tencellub.head$Peptide.Name), ]
tencellub.head[grep("HLA-B ubi K340", tencellub.head$Peptide.Name)[1], "Peptide.Name"] <- "HLA-B.7 ubi Y344"
tencellub.head[grep("HLA-B ubi K340", tencellub.head$Peptide.Name)[1], "Peptide.Name"] <- "HLA-B.4 ubi Y344"
tencellub.head[grep("HLA-C", tencellub.head$Peptide.Name), ]
tencellub.head[grep("HLA-C ubi K365", tencellub.head$Peptide.Name)[1], "Peptide.Name"] <- "HLA-C.4 ubi  K365"
any(duplicated(tencellub.head$Peptide.Name)) #  FALSE!
# Okay, proceed
tencellubdata <- tencellub[,grep("Intensity", names(tencellub))]
dim(tencellubdata) # 4683   26
names(tencellubdata) <- sapply(names(tencellubdata), function (x) unlist(strsplit(x, "Intensity."))[2])

# make zero into NA, which it is.
zer0.tc <- which(tencellubdata==0, arr.ind = TRUE)
tencellubdata <- replace (tencellubdata, zer0.tc, NA)
nmissing(tencellubdata)/(dim(tencellubdata)[1]*dim(tencellubdata)[2])
# 78% NA
# To make column names identical, remove dates on the left and the last three characters
names(tencellubdata) <- sapply(names(tencellubdata), function (x) unlist(strsplit(x, "GZ_"))[2])
tcdnames <- names(tencellubdata)
tcdnames <- make.unique(str_sub(tcdnames, 1, str_length(tcdnames)-3))
names(tencellubdata) <- tcdnames
identical(names(tencellphosdata), names(tencellubdata))  # FALSE
data.frame(names(tencellphosdata), names(tencellackdata), names(tencellubdata))
# Column names are out of order!
# add rownames to data
rownames(tencellubdata) <- tencellub.head$Peptide.Name
# Assume columns are out of order, not just the column names
names(tencellphosdata)[c(1:14,16,15,17:26)] # check
tencellphosdata.1 <- tencellphosdata[, c(1:14,16,15,17:26)]
data.frame(names(tencellphosdata.1), names(tencellackdata), names(tencellubdata))
# Fix typos
names(tencellphosdata.1)[15] <- "H2286_Dasatinib"
names(tencellackdata)[19] <- "H366_Dasatinib"
names(tencellackdata)[20] <- "H366_DMSO"
tencelldata <- rbind(tencellphosdata.1, tencellackdata, tencellubdata)
dim(tencelldata)
# [1] 10624    26
any(is.na(rownames(tencelldata))) # F
# How many PTMs are detected more than twice?
tencellsamples <- apply(tencelldata, 1, filled)
hist(tencellsamples, breaks=50, col="magenta")
tencelltrimmed <- tencelldata[which(tencellsamples>2),]
dim(tencelltrimmed) # now 7847
#######################################################################################
#save(tencellack, tencellack.head, tencellackdata, tencellackname, tencellpath, tencellphos, tencellphos.head, tencellphosdata, tencellphosdata.1, tencellphosname, tencellub, tencellub.head, tencellubname, tencelldata, file=paste("_LINCS/_KarenGuolin/", "TenCell.RData", sep=""))
#######################################################################################
# Ratios
# What is the overlap between data colums used to calculate ratios?
#
H2228CrizRatio <- tencelldata[, "H2228_Crizotinib"]/tencelldata[, "H2228_DMSO"]
nmissing (tencelldata[, "H2228_Crizotinib"]); nmissing(tencelldata[, "H2228_DMSO"]); nmissing (H2228CrizRatio)
filled (tencelldata[, "H2228_Crizotinib"]); filled(tencelldata[, "H2228_DMSO"]); filled (H2228CrizRatio)
filled (tencelldata[, "H2228_Crizotinib"])/length(tencelldata[, "H2228_Crizotinib"]); filled(tencelldata[, "H2228_DMSO"])/length(tencelldata[, "H2228_DMSO"]); filled (H2228CrizRatio)/length(H2228CrizRatio)
# There is a loss of about a third of the data
H2228CrizotinibRatio <- tencelldata[, "H2228_Crizotinib"]/tencelldata[, "H2228_DMSO"]
H3122CrizotinibRatio <- tencelldata[, "H3122_Crizotinib"]/tencelldata[, "H3122_DMSO"]
H2228CrizotinibRatio <- tencelldata[, "H2228_Crizotinib"]/tencelldata[, "H2228_DMSO"]
HCC4006_ErlotinibRatio <- tencelldata[, "HCC4006_Erlotinib"]/tencelldata[, "HCC4006_DMSO"]
HCC78_CrizotinibRatio <- tencelldata[, "HCC78_Crizotinib"]/tencelldata[, "HCC78_DMSO"]
HCC827_ErlotinibRatio <- tencelldata[, "HCC827_Erlotinib"]/tencelldata[, "HCC827_DMSO"]
PC9_ErlotinibRatio <- tencelldata[, "PC9_Erlotinib"]/tencelldata[, "PC9_DMSO"]
H1781_AfatinibRatio <- tencelldata[, "H1781_Afatinib"]/tencelldata[, "H1781_DMSO"]
H2286_DasatinibRatio <- tencelldata[, "H2286_Dasatinib"]/tencelldata[, "H2286_DMSO"]
H3122_Crizotinib.1Ratio <- tencelldata[, "H3122_Crizotinib.1"]/tencelldata[, "H3122_DMSO.1"]
H366_DasatinibRatio <- tencelldata[, "H366_Dasatinib"]/tencelldata[, "H366_DMSO"]
HCC78_Crizotinib.1Ratio <- tencelldata[, "HCC78_Crizotinib.1"]/tencelldata[, "HCC78_DMSO.1"]
PC9_Erlotinib.1Ratio <- tencelldata[, "PC9_Erlotinib.1"]/tencelldata[, "PC9_DMSO.1"]
STE.1_CrizotinibRatio <- tencelldata[, "STE.1_Crizotinib"]/tencelldata[, "STE.1_DMSO"]
tencellratios <- data.frame(H1781_AfatinibRatio, H2228CrizotinibRatio, H2286_DasatinibRatio, H3122_Crizotinib.1Ratio, H3122CrizotinibRatio, H366_DasatinibRatio, HCC4006_ErlotinibRatio, HCC78_Crizotinib.1Ratio, HCC78_CrizotinibRatio, HCC827_ErlotinibRatio, PC9_Erlotinib.1Ratio, PC9_ErlotinibRatio, STE.1_CrizotinibRatio)
#
rownames(tencellratios) <- rownames(tencelldata)
identical(rownames(tencelldata), rownames(tencellratios))
# TRUE
# Make limits as for ld data
hi.ratio <- which(tencellratios>=100, arr.ind = TRUE)
low.ratio <- which(tencellratios <=1/100, arr.ind = TRUE)
tencellratios.lim <- replace (tencellratios, hi.ratio, 100) 
tencellratios.lim <- replace (tencellratios.lim, low.ratio, 1/100) 
# log2
tencellratios.log2 <- log2(tencellratios)
tencellratios.lim.log2 <- log2(tencellratios.lim)
nmissing(tencellratios.lim.log2)/(dim(tencellratios.lim.log2)[1]*dim(tencellratios.lim.log2)[2])
# [1] 84% NA 
identical(rownames(tencellratios.lim.log2), rownames(tencelldata))
# [1] TRUE
boxplot(tencellratios.lim.log2)
#
tencelldata.log2 <- log2(tencelldata)
boxplot(tencelldata.log2)
nmissing(tencelldata.log2)/(dim(tencelldata.log2)[1]*dim(tencelldata.log2)[2])
# 74% NA - about 10% data is lost from ratios
identical(rownames(tencellratios.lim.log2), rownames(tencelldata.log2))
# True
# Check for NA in names
whichrowNA <- which(grepl("NA ", rownames(tencellratios.lim.log2)))
rownames(tencellratios.lim.log2)[whichrowNA] 
# "NA ubi K162; NA ubi K132" "NA ubi K359" "NA ubi K252; NA ubi K338; ...  
whichrowNAd <- which(grepl("NA ", rownames(tencelldata.log2)))
rownames(tencelldata.log2)[whichrowNAd]    
# NA ubi K162; NA ubi K132"    "NA ubi K252; NA ubi K338;...[66] "NA ubi K359"    
badnames <- c("NA ubi K162; NA ubi K132", "NA ubi K359", "NA ubi K252; NA ubi K338; TUBB4A ubi K252; TUBB ubi K252; TUBB4B ubi K252; TUBB3 ubi K252; TUBB2A ubi K252; TUBB8 ubi K252; TUBB6 ubi K252; TUBB2B ubi K252")
    
# Make trimmed version
tencellratios.lim.log2.trimmed <- tencellratios.lim.log2[which(tencellsamples>2),]
#######################################################################################
# Merge with gzdata <- cbind(kglog2data, kgratios.lim.log2)
gztencelldata <- cbind(tencelldata.log2, tencellratios.lim.log2)
gzdata.all <- merge(gzdata, gztencelldata, by="row.names", all=TRUE)
dim(gzdata.all)
rownames(gzdata.all) <- gzdata.all$Row.names
gzdata.allz <- gzdata.all[,2:ncol(gzdata.all)]
nmissing(gzdata.allz)/(dim(gzdata.allz)[1]*dim(gzdata.allz)[2])
# 83% NA
gztencelldata.trimmed <- cbind(log2(tencelltrimmed), tencellratios.lim.log2.trimmed)
gzdata.all.trimmed <- merge(gzdata, gztencelldata.trimmed, by="row.names", all=TRUE)
whichrowNA <- which(grepl("NA ", rownames(gzdata.allt))
                    rownames(gzdata.allt)[whichrowNA]                    
                    
rownames(gzdata.all.trimmed) <- gzdata.all.trimmed$Row.names
gzdata.allt <- gzdata.all.trimmed[,2:ncol(gzdata.all.trimmed)]
nmissing(gzdata.allt)/(dim(gzdata.allt)[1]*dim(gzdata.allt)[2])
# now 78% NA; tencell trimmed is 71%
#######################################################################################
# save(tencellack, tencellack.head, tencellackdata, tencellackname, tencellpath, tencellphos, tencellphos.head, tencellphosdata, tencellphosdata.1, tencellphosname, tencellub, tencellub.head, tencellubname, tencelldata, tencelldata.log2, tencellratios.lim, tencellratios.log2, tencellratios.lim.log2, gzdata.all, gzdata.allz, file=paste("_LINCS/_KarenGuolin/", "TenCell.RData", sep=""))
#######################################################################################
# Calculate dissimilarity matrices on gz all data data next: 
# using gzdata.allz (gzall) data 
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Do first with tencell trimmed, then rename: 
#gztencelltrimmed.cor <- gzallt.cor
#dissimilarity.gztencell <- dissimilarity.gzallt
#diss.gztencell.noabs <- diss.gzallt.noabs
#gztencell.dist <- gzallt.dist
#gzatencell.sed <- gzallt.sed
#eu.gztencell.tsne <- eu.gzallt.tsne
#sp.gztencell.tsne <- sp.gzallt.tsne
#sed.gztencell.tsne <- sed.gzallt.tsne
#eu.sp.sed.gztencell <- eu.sp.sed.gzallt
#eu.sp.sed.gztencell.data <- eu.sp.sed.gzallt.data
#XXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Spearman
gzall.cor <- cor(t(gzdata.allz), use = "pairwise.complete.obs", method = "spearman")
# The trimmed data is throwing errors "Error: vector memory exhausted (limit reached?)" on office computer
# so try Pearson
# Try also on laptop. Allt = gzdata combined with tencell trimmed. Works!
gzallt.cor <- cor(t(gzdata.allt), use = "pairwise.complete.obs", method = "spearman")
# Pearson Works, though there were many zero std. deviations
diag(gzall.cor) <- NA
diag(gzallt.cor) <- NA
dissimilarity.gzall <- 1 - abs(gzall.cor)
diss.gzall.noabs <- 1 - gzall.cor
dissimilarity.gzallt <- 1 - abs(gzallt.cor)
diss.gzallt.noabs <- 1 - gzallt.cor
# set NA to two orders of magnitude higher than max distance
dissimilarity.gzall[is.na(dissimilarity.gzall)] <- 100*max(dissimilarity.gzall, na.rm=T)
diss.gzall.noabs[is.na(diss.gzall.noabs)] <- 50*max(diss.gzall.noabs, na.rm=T) 
dissimilarity.gzallt[is.na(dissimilarity.gzallt)] <- 100*max(dissimilarity.gzallt, na.rm=T)
diss.gzallt.noabs[is.na(diss.gzallt.noabs)] <- 50*max(diss.gzallt.noabs, na.rm=T) 
# check
max(dissimilarity.gzall)
max(diss.gzall.noabs)
max(dissimilarity.gzallt)
max(diss.gzallt.noabs)
# now max=100
# Euclid	
gzall.dist = as.matrix (dist (gzdata.allz), method = "euclidean")  
gzallt.dist = as.matrix (dist (gzdata.allt), method = "euclidean")  
# check 
max.na(gzall.dist)	# Check for Inf values here
max.na(gzallt.dist)	# Check for Inf values here
gzall.dist[is.na(gzall.dist)] <- 100*max(gzall.dist, na.rm=T)
gzall.dist.1 <- 100*gzall.dist/max(gzall.dist, na.rm=T)  # now max=100
gzallt.dist[is.na(gzallt.dist)] <- 100*max(gzallt.dist, na.rm=T)
gzallt.dist.1 <- 100*gzallt.dist/max(gzallt.dist, na.rm=T)  # now max=100
# SED: combinde Euclid and Spearman w/o taking absolute value.	
gzall.sed <- (gzall.dist.1 + diss.gzall.noabs)/2
gzallt.sed <- (gzallt.dist.1 + diss.gzallt.noabs)/2
#
# Save these because they took a long time: 2.17 GB
#save(gzall.cor, gzallt.cor, dissimilarity.gzall, diss.gzall.noabs, gzall.dist, gzall.dist.1, gzall.sed, file="/Users/mark_grimes/Software\ Images/GZTenCellMatrices.RData")
# >>>>>-------> Fix to include only PTMs from >2 experiments for CCCN construction
if(dim(gzallt.cor)[1]==9400) {
  gzallt.cor.1 <- gzallt.cor[rownames(gzallt.cor) %in% rownames(gzdata.allt), colnames(gzallt.cor) %in% rownames(gzdata.allt)]  }
gzallt.cor <- gzallt.cor.1
# List elements taken care of below

##___________________________________
#
# Rtsne: Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding
require(Rtsne)
#
# Do with perplexity = 15 for slightly more resolution; theta = 0.25 for slightly more accuracy, max_iter=5000, for stabiliztion of groups
eu.gzall.tsne.list <- Rtsne(as.matrix(gzall.dist), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca=FALSE)
eu.gzallt.tsne.list <- Rtsne(as.matrix(gzallt.dist), dims = 3, perplexity = 15, theta = 0.25, max_iter=5000, check_duplicates = FALSE, pca=FALSE)
#    
eu.gzall.tsne <- eu.gzall.tsne.list$Y
eu.gzallt.tsne <- eu.gzallt.tsne.list$Y
#
sp.gzall.tsne.list <- Rtsne(as.matrix(dissimilarity.gzall), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca=FALSE)
sp.gzallt.tsne.list <- Rtsne(as.matrix(dissimilarity.gzallt), dims = 3, perplexity = 15, theta = 0.25, max_iter=5000, check_duplicates = FALSE, pca=FALSE)
#   
sp.gzall.tsne <- sp.gzall.tsne.list$Y
sp.gzallt.tsne <- sp.gzallt.tsne.list$Y
# SED without abs(cor)
sed.gzall.tsne.list <- Rtsne(as.matrix(gzall.sed), dims = 3, perplexity = 15, theta = 0.25, max_iter=5000, check_duplicates = FALSE, pca=FALSE)
sed.gzallt.tsne.list <- Rtsne(as.matrix(gzallt.sed), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca=FALSE)
#   
sed.gzall.tsne <- sed.gzall.tsne.list$Y
sed.gzallt.tsne <- sed.gzallt.tsne.list$Y
#    
#######################################################################################
# save(tencellack, tencellack.head, tencellackdata, tencellackname, tencellpath, tencellphos, tencellphos.head, tencellphosdata, tencellphosdata.1, tencellphosname, tencellub, tencellub.head, tencellubname, tencelldata, tencelldata.log2, tencellratios.lim, tencellratios.log2, tencellratios.lim.log2, gzdata.all, gzdata.allz, gzdata.all.raw, eu.gzall.tsne, sp.gzall.tsne, sed.gzall.tsne, gzdata.all.trimmed, gztencelldata.trimmed, tencellratios.lim.log2.trimmed, tencelltrimmed, gzdata.allt, eu.gzallt.tsne, sp.gzallt.tsne, sed.gzallt.tsne, file=paste("_LINCS/_KarenGuolin/", "TenCell.RData", sep=""))
#######################################################################################
#
# Without ratios
#
gzdata.all.raw <- merge(kglog2data, tencelldata.log2, by="row.names", all=TRUE)
rownames(gzdata.all.raw) <- gzdata.all.raw$Row.names
gzdata.all.raw <- gzdata.all.raw[,2:ncol(gzdata.all.raw)]
#
#_________________________________________________
# Try without ratio data
#_________________________________________________
#######################################################################################
# Calculate dissimilarity matrices on gz all data data next: 
# using gzdata.all.raw (gzdata.all.raw) data 
# NOTE: Many errors and warnings here! - Don't use
# Spearman
gzall.raw.cor <- cor(t(gzdata.all.raw), use = "pairwise.complete.obs", method = "spearman")
diag(gzall.raw.cor) <- NA
hist(unlist(gzall.raw.cor), breaks=100, col="blue")
dissimilarity.gzall.raw <- 1 - abs(gzall.raw.cor)
diss.gzall.raw.noabs <- 1 - gzall.raw.cor
# set NA to two orders of magnitude higher than max distance
dissimilarity.gzall.raw[is.na(dissimilarity.gzall.raw)] <- 100*max(dissimilarity.gzall.raw, na.rm=T)
diss.gzall.raw.noabs[is.na(diss.gzall.raw.noabs)] <- 50*max(diss.gzall.raw.noabs, na.rm=T) 
# check
max(dissimilarity.gzall.raw)
max(diss.gzall.raw.noabs)
# now max=100
# Euclid	
gzall.raw.dist = as.matrix (dist (gzdata.all.raw), method = "euclidean")  
# check 
max.na(gzall.raw.dist)	# Check for Inf values here
hist(unlist(gzall.raw.dist), breaks=100, col="yellow")
#
gzall.raw.dist[is.na(gzall.raw.dist)] <- 100*max(gzall.raw.dist, na.rm=T)
gzall.raw.dist.1 <- 100*gzall.raw.dist/max(gzall.raw.dist, na.rm=T)  # now max=100
# SED: combinde Euclid and Spearman w/o taking absolute value.	
gzall.raw.sed <- (gzall.raw.dist.1 + diss.gzall.raw.noabs)/2
#
# Save these because they took a long time:  GB
#save(gzall.cor, dissimilarity.gzall, diss.gzall.noabs, gzall.dist, gzall.dist.1, gzall.sed, gzall.raw.cor, dissimilarity.gzall.raw, diss.gzall.raw.noabs, gzall.raw.dist, gzall.raw.dist.1, gzall.raw.sed, file="/Users/mark_grimes/Software\ Images/GZTenCellMatrices.RData")
#save(gzall.cor, dissimilarity.gzall, diss.gzall.noabs, gzall.dist, gzall.dist.1, gzall.sed, gzall.raw.cor, dissimilarity.gzall.raw, diss.gzall.raw.noabs, gzall.raw.dist, gzall.raw.dist.1, gzall.raw.sed, file="/Volumes/Terra_Byte/R_Archive_2/_LINCS/GZTenCellMatrices.RData")
##___________________________________
#
# Rtsne: Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding
require(Rtsne)
#
# Do with perplexity = 15 for slightly more resolution; theta = 0.25 for slightly more accuracy
eu.gzall.raw.tsne.list <- Rtsne(as.matrix(gzall.raw.dist), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca=FALSE)
#    
eu.gzall.raw.tsne <- eu.gzall.raw.tsne.list$Y
#
sp.gzall.raw.tsne.list <- Rtsne(as.matrix(dissimilarity.gzall.raw), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca=FALSE)
#   Error: vector memory exhausted (limit reached?). Try with PCA
sp.gzall.raw.tsne.list <- Rtsne(as.matrix(dissimilarity.gzall.raw), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca=TRUE)
# Same error (?)

sp.gzall.raw.tsne <- sp.gzall.raw.tsne.list$Y
# SED without abs(cor)
sed.gzall.raw.tsne.list <- Rtsne(as.matrix(gzall.raw.sed), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca=FALSE)
#  Error: vector memory exhausted (limit reached?) 
sed.gzall.raw.tsne <- sed.gzall.raw.tsne.list$Y
#    t-SNE fails without including ratios!
# save(gzall.cor, gztencelltrimmed.cor, gzallt.cor, dissimilarity.gzall, dissimilarity.gztencell, diss.gzall.noabs, diss.gztencell.noabs, gzall.dist, gzall.dist.1, gztencell.dist, gzall.sed, gzatencell.sed, gzall.raw.cor, dissimilarity.gzall.raw, diss.gzall.raw.noabs, gzall.raw.dist, gzall.raw.dist.1, gzall.raw.sed, eu.gzall.tsne, sp.gzall.tsne, sed.gzall.tsne, file="/Volumes/Terra_Byte/R_Archive_2/_LINCS/GZTenCellMatrices2.RData")
save(gzall.cor, gztencelltrimmed.cor, gzallt.dist, gzallt.cor, dissimilarity.gzall, dissimilarity.gzallt, diss.gzallt.noabs, dissimilarity.gztencell, diss.gzall.noabs, diss.gztencell.noabs, gzall.dist, gzall.dist.1, gztencell.dist, gzall.sed, gzatencell.sed, gzallt.dist.1, gzallt.sed, gzall.raw.cor, dissimilarity.gzall.raw, diss.gzall.raw.noabs, gzall.raw.dist, gzall.raw.dist.1, gzall.raw.sed, eu.gzall.tsne, sp.gzall.tsne, sed.gzall.tsne, eu.gzallt.tsne, sp.gzallt.tsne, sed.gzallt.tsne, file="/Volumes/Terra_Byte/R_Archive_2/_LINCS/GZTenCellMatrices3.RData")

# save(gzall.cor, gzallt.cor, gzallt.dist, dissimilarity.gzall, dissimilarity.gzallt, diss.gzall.noabs, gzall.dist, gzall.dist.1, gzall.sed, diss.gzallt.noabs, gzallt.dist, gzallt.dist.1, gzallt.sed, gzall.raw.cor, dissimilarity.gzall.raw, diss.gzall.raw.noabs, gzall.raw.dist, gzall.raw.dist.1, gzall.raw.sed, eu.gzall.tsne, sp.gzall.tsne, sed.gzall.tsne, file="/Users/mark_grimes/Software\ Images/GZTenCellMatrices.RData")

# Three dimensional plots
require(rgl)
plot3d(eu.gzall.tsne, type="s", radius=0.76, col="forestgreen")
plot3d(sp.gzall.raw.tsne, type="s", radius=0.76, col="green")
plot3d(sed.gzall.tsne, , type="s", radius=0.76, col="red")
# New
plot3d(eu.gzallt.tsne, type="s", radius=1.5, col="forestgreen")
plot3d(sp.gzallt.tsne, type="s", radius=1.5, col="green")
plot3d(sed.gzallt.tsne, , type="s", radius=1.5, col="red")

#######################################################################################
eu.gzall.list <- make.clusterlist(eu.gzall.tsne, 3.5, gzdata.allz) # 216 groups
esizes.gzall <- sapply(eu.gzall.list, function(x) dim(x)[1])
sp.gzall.list <- make.clusterlist(sp.gzall.tsne, 3.5, gzdata.allz) # 171 groups
spsizes.gzall <- sapply(sp.gzall.list, function(x) dim(x)[1])
sed.gzall.list <- make.clusterlist(sed.gzall.tsne, 3.5, gzdata.allz) # 224 groups
sedsizes.gzall <- sapply(sed.gzall.list, function(x) dim(x)[1])
#
eu.gzallt.list <- make.clusterlist(eu.gzallt.tsne, 3.8, gzdata.allt) 
# 261 groups round 1; 497 round 2 @ 3.5; 415 @ 3.8
esizes.gzallt <- sapply(eu.gzallt.list, function(x) dim(x)[1])
sp.gzallt.list <- make.clusterlist(sp.gzallt.tsne, 3.8, gzdata.allt) 
# 290 groups round 1;  593 round 2 @ 3.5; 504 @ 3.8
spsizes.gzallt <- sapply(sp.gzallt.list, function(x) dim(x)[1])
sed.gzallt.list <- make.clusterlist(sed.gzallt.tsne, 3.0, gzdata.allt) 
# 284 groups round 1: 266 round 2 @ 3.5; 345 @ 3.0
sedsizes.gzallt <- sapply(sed.gzallt.list, function(x) dim(x)[1])
#
hist(esizes.gzall, breaks=100, col="green") 
# note several large clusters that should be broken up
hist(spsizes.gzall, breaks=100, col="yellow") 
# note several large clusters that should be broken up
hist(sedsizes.gzall, breaks=100, col="blue") 
# note several large clusters that should be broken up
hist(esizes.gzallt, breaks=100, col="red")
hist(spsizes.gzallt, breaks=100, col="blue")
hist(sedsizes.gzallt, breaks=100, col="green")

# Examine the large group data
fractNA <- function(df) {
  result <- nmissing(df)/(dim(df)[1]*dim(df)[2])
  return(result)
}


###########----------
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
##
# 
#  Concatenate group list files
#		make group names unique
eu.gzallt.df <- ldply(eu.gzallt.list)[,2:3]
sp.gzallt.df <- ldply(sp.gzallt.list)[,2:3]
sed.gzallt.df <- ldply(sed.gzallt.list)[,2:3]	# Further partition large groups
eu.gzallt.df $group <- paste(noquote(eu.gzallt.df $group), noquote("e"), sep="", collapse=NULL)
sp.gzallt.df $group <- paste(noquote(sp.gzallt.df $group), noquote("s"), sep="", collapse=NULL)
sed.gzallt.df $group <- paste(noquote(sed.gzallt.df $group), noquote("sed"), sep="", collapse=NULL)
gzalltgroups.df <- rbind(eu.gzallt.df, sed.gzallt.df, sp.gzallt.df)
# 41394 round 1; 28200 round 2
gzalltgroups.tab <- table(gzalltgroups.df)	
#  dim 7847 peptides X 835 groups; round 2 9400 x 1435 (1264 round 3 with tweaked toolong)
#
eu.gzallt.genes <- lapply(eu.gzallt.list, extract.genes.from.clist)
sp.gzallt.genes <- lapply(sp.gzallt.list, extract.genes.from.clist)
sed.gzallt.genes <- lapply(sed.gzallt.list, extract.genes.from.clist)
#
eu.gzallt.peps <- lapply(eu.gzallt.list, extract.peps.from.clist)
sp.gzallt.peps <- lapply(sp.gzallt.list, extract.peps.from.clist)
sed.gzallt.peps <- lapply(sed.gzallt.list, extract.peps.from.clist)
#
eu.sp.gzallt <- list.common(eu.gzallt.peps, sp.gzallt.peps, keeplength=2)
eu.sp.gzallt.sizes <- sapply(eu.sp.gzallt, length)
eu.sp.sed.gzallt <- list.common(eu.sp.gzallt, sed.gzallt.peps, keeplength=2)
eu.sp.sed.gzallt.sizes <- sapply(eu.sp.sed.gzallt, length) 
hist(eu.sp.sed.gzallt.sizes, breaks=100, col="gold")
# max = 111; length = 612; round 2 max 88, length 878; 111 and 839 round 3
# patch to deal with untrimmed data
clust.data.from.vec <- function(vec, tbl) {
  if(class(vec)=="list") {vec <- unlist(vec)}
  at <- tbl[vec,]
  acol <- names(at[,which(numcolwise(filled)(at) != 0)])
  if(length(acol)  == 1) {
    ats <- data.frame(cbind (rownames(at), as.numeric(at[, acol])))
    names(ats) <- c("Gene.Name", acol)
  }
  if(length(acol) >= 2) {
    ats <- cbind(rownames(at), at[, acol])
    names(ats)[1] <- "Gene.Name" } 
  clust.data <- ats
  return (clust.data)	}
# # >>>>>
# Examine all clusters
eu.sp.sed.gzallt.data <- lapply(eu.sp.sed.gzallt, clust.data.from.vec, tbl= gzdata.allt) 
# This works the second time around. First time produces errors. try this:
# eu.sp.sed.gzallt.data <- list()
# for (i in 1:length(eu.sp.sed.gzallt)) {
  if (length(intersect(eu.sp.sed.gzallt[[i]], rownames(gzdata.allt)))==0) next
  at <- gzdata.allt[unlist(eu.sp.sed.gzallt[[i]]),]
  if(dim(at)[1]<2 | dim(at)[2]<2) next
  eu.sp.sed.gzallt.data[[i]] <- clust.data.from.vec(eu.sp.sed.gzallt[[i]], tbl=gzdata.allt)
  print(i)
}
# error at 243 (fixed above)
# Note : names missing.
gzdata.allt[rownames(gzdata.allt) %in% eu.sp.sed.gzallt[[i]],]
# data in only one column! xxxx
# first data set was not trimmed for PTMs are detected more than twice
alltsamples <- apply(gzdata.allt, 1, filled)
hist(alltsamples, breaks=50, col="magenta")
gzdata.allt.t <- gzdata.allt[which(alltsamples>2),]
dim(gzdata.allt.t) # now 9215 from 9400
gzdata.allt <- gzdata.allt.t
# Go to above 
# Print out and delete bad ones
bad.clusterlist <- list()
for (i in 1:length(eu.sp.sed.gzallt)) {
  if (length(intersect(eu.sp.sed.gzallt[[i]], rownames(gzdata.allt)))==0)  {
    print(i)
      bad.clusterlist[[i]] }
}
length(eu.sp.sed.gzallt.data) # 839

# Find any that may be included in clusters
badptms <- unique(outersect(rownames(gzdata.all), rownames(gzdata.allt))) # 2844
testnames <- unique(unlist(eu.sp.sed.gzallt))
length(intersect (badptms, testnames))  # 181
essgzallt.1 <- lapply(eu.sp.sed.gzallt, function(x) x %w/o% badptms)
essgzallt.1.sizes <- sapply(essgzallt.1, length)
removed <- eu.sp.sed.gzallt.sizes - essgzallt.1.sizes
removed <- removed[removed>0]
# note: some entire clusters removed
essgzallt.2 <- essgzallt.1[essgzallt.1.sizes>0] # now 818 clusters from 839
essgzallt <- essgzallt.2
# >>>> Now redo data list
essgzallt.data <- lapply(essgzallt, clust.data.from.vec, tbl= gzdata.allt) 

#

# _________________________
# # # >>>>>
# Evalauate clusters
# use lincsclust.eval; clusterlist=essgzallt.data; tbl.sc=gzdata.allt
lincsclust.eval <- function(clusterlist, tbl.sc) {
  evaluation <- data.frame(0)
  names(evaluation)[1] <- "Group"
  key  <- data.frame(1:length(rownames(tbl.sc)))
  key$ptm.Name <- rownames(tbl.sc)
  for (i in 1:length(clusterlist)) {
    cat("Starting Group", i, "\n")
    evaluation[i,1] <- i
    evaluation$Group.Name[i] <- names(clusterlist)[i]
    #
    evaluation$no.ptms[i] <- length(clusterlist[[i]]$ptm.Name)
    if(length(clusterlist[[i]]$ptm.Name) == 1) { 
      at = data.frame(tbl.sc[clusterlist[[i]]$ptm.Name, ])
    } else {
      at = data.frame(tbl.sc[key$ptm.Name %in% clusterlist[[i]]$ptm.Name, ]) }
    # get rid of ratios for evaluation calculations and take absolute value
    if(any (grepl("atio", names(at)))) at = abs(at [,-grep("atio", names(at))])
    # previous use: at <- at[-which(apply(at, 1, filled) == 0),]
    # better: at[, which(numcolwise(filled)(at) != 0)]
    # names() doesn't work with single column
    if (length(which(numcolwise(filled)(at) != 0)) > 1) {
      acol <- names(at[,which(numcolwise(filled)(at) != 0)])  
      evaluation$no.samples[i] <- length(acol)
      at <- at[, acol] } else evaluation$no.samples[i] <- 1
    evaluation$total.signal[i] <- sum(abs(at), na.rm=TRUE)
    if (length(which(numcolwise(filled)(at) != 0)) == 1 || length(clusterlist[[i]]$ptm.Name) == 1) {
      evaluation$culled.by.slope[i] <- length(clusterlist[[i]]$ptm.Name) 
      evaluation$percent.NA[i] <- 0
      evaluation$percent.singlesampleptms[i] <- 100
      evaluation$percent.singleptmsamples[i] <- 100
    } else	{
      evaluation$percent.NA[i] <-  100*(sum(numcolwise(nmissing)(at)) / (dim(at)[1]*dim(at)[2]))
      #singlesampleptms <- at[, which(numcolwise(filled)(at) == 1 )]
      singlesampleptms <- at[which(apply(at, 1, filled) == 1),]
      evaluation$percent.singlesampleptms[i] <- 100*(nrow(singlesampleptms) / dim(at)[1]) 
      singleptmsamples <- sum(numcolwise(filled)(at) == 1)
      evaluation$percent.singleptmsamples[i] <- 100*(singleptmsamples/dim(at)[2])
      cluster.mo <- at[order(-as.vector(colwise(sum.na)(data.frame(t(abs(at)))))), order(-as.vector(numcolwise(sum.na)(data.frame(abs(at)))))]
      slope <- apply(cluster.mo, 1, get.slope.a)
      badslope <- c(names(which(is.na(slope))), names(which(slope > 0)))
      evaluation$culled.by.slope[i] <- length(badslope)
      #
      cat("\n", length(badslope), "ptms culled by slope", "\n")
    }		}	 
  #  Total signal scaled to percent NA = intensity
  clearptms <- evaluation$no.ptms - evaluation$culled.by.slope # may be 0
  realsamples <- evaluation$no.samples - (evaluation$no.samples * evaluation$percent.singleptmsamples/100) # may be 0
  intensity <- evaluation$total.signal - (evaluation$total.signal * evaluation$percent.NA/100)
  # calibrate intensity according to real samples and clear ptms
  # - goal is to reward a high density of appropriate data
  evaluation$intensity <- intensity
  evaluation$Index  <- ((1 + realsamples) * (1 + clearptms) / (1 + evaluation$percent.NA))/evaluation$no.ptms
  eval.sort <- evaluation[order(-evaluation$Index, evaluation$percent.NA), c("Group", "Group.Name", "no.ptms",  "culled.by.slope", "percent.singlesampleptms","no.samples", "percent.singleptmsamples", "total.signal", "percent.NA", "intensity", "Index" )] 
  return(eval.sort)	
}
gzclust.eval.df.1 <- lincsclust.eval(eu.sp.sed.gzallt.data, tbl.sc=gzdata.allt)
# 143 is bad (NA) ;242 is bad; 426, 431, 536, 611, 640, 644...manually remove & redo above after deleting bad ones
# OR
# Better: 
# Try with already trimmed data above
# essgzallt.data <- lapply(essgzallt, clust.data.from.vec, tbl= gzdata.allt) 

gzclust.eval.df.test <- lincsclust.eval(essgzallt.data, tbl.sc=gzdata.allt)
# Test new version with names changed
gzclust.eval.df <- ptmsclust.eval(essgzallt.data, tbl.sc=gzdata.allt)
# Warnings but no NA clusters
# loop to examine heatmaps
for (i in 818:810){
  graph.clust6d.l(essgzallt.data[[gzclust.eval.df$Group.Name[i]]])
}
graph.clust6d.l(essgzallt.data[[gzclust.eval.df$Group.Name[818]]])
# As noted above: 
# In about a third of the data from either treatment or control PTMs missing.
# Find clusters that have two or more ratio columns.
essgzallt.data.ratios <- lapply(essgzallt.data, function (x) x[, grepl("atio", names(x)), drop=FALSE])
edr.cols <- sapply(essgzallt.data.ratios, function (x) dim(x)[2])
hist(unlist(edr.cols), breaks=20, col="turquoise")
range(edr.cols)
essgzallt.data.ratios <-  essgzallt.data.ratios[which(unlist(edr.cols)>0)] # 767


gzclust.eval.df$contains.ratio.data <- sapply(gzclust.eval.df$Group.Name, function (y) y %in% names(edr))
# remove this and just add no. ratio cols.
gzclust.eval.df <- gzclust.eval.df[ , names(gzclust.eval.df) %w/o% "contains.ratio.data"]
gzclust.eval.df$no.ratio.cols <- edr.cols
write.table(gzclust.eval.df, file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "gzclust.eval.df.txt", sep=""), row.names = FALSE, sep="\t")

# Evaluate ratio clusters with more than three data columns with modified function
essgzallt.data.ratios.3cols <- essgzallt.data.ratios[which(unlist(edr.cols)>2)] # 524
ratioclust.eval.df <- ratioclust.eval(clusterlist=essgzallt.data.ratios.3cols)
write.table(ratioclust.eval.df, file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "gzratioclusters.eval.df.txt", sep=""),row.names = FALSE, sep="\t")
graph.clust6d.l(edr["118.227.333"])
# loop to examine heatmaps
for (i in 1:10){
  graph.clust6d.l(edr[[ratioclust.eval.df$Group.Name[i]]])
}
# _
graph.clust6d.l(edr[["27.282.215"]])
graph.clust6d.l(essgzallt.data[["27.282.215"]])
graph.clust6d.l(essgzallt.data.ratios[["27.282.215"]])
graph.clust6d.l(essgzallt.data.ratios[["2.58.53"]])
graph.clust6d.l(essgzallt.data[["1.225.171"]])

#________________________

save(tencellack, tencellack.head, tencellackdata, tencellackname, tencellpath, tencellphos, tencellphos.head, tencellphosdata, tencellphosdata.1, tencellphosname, tencellub, tencellub.head, tencellubname, tencelldata, tencelldata.log2, tencellratios.lim, tencellratios.log2, tencellratios.lim.log2, gzdata.all, gzdata.allz, gzdata.all.raw, eu.gzall.tsne, sp.gzall.tsne, sed.gzall.tsne, gzdata.all.trimmed, gztencelldata.trimmed, tencellratios.lim.log2.trimmed, tencelltrimmed, gzdata.allt, eu.gzallt.tsne, sp.gzallt.tsne, sed.gzallt.tsne, eu.gzall.list, sp.gzall.list, sed.gzall.list, esizes.gzall, spsizes.gzall, sedsizes.gzall, eu.gzallt.list, sp.gzallt.list, sed.gzallt.list, esizes.gzallt, eu.sp.sed.gzallt, eu.sp.sed.gzallt.sizes, eu.sp.sed.gzallt.data, spsizes.gzallt, eu.gztencell.tsne, sp.gztencell.tsne, sed.gztencell.tsne, sedsizes.gzallt, eu.sp.sed.gztencell, eu.sp.sed.gztencell.data, essgzallt, essgzallt.data, essgzallt.data.ratios, gzclust.eval.df, ratioclust.eval.df, gzalltgenecccn.edges, gzallt.gene.cccn.g, gzallt.gene.cccn0, gzallt.gene.cccn.na, gzallt.cccn.g, gzallt.cccn, pepcorredges.dual.neg.v1, pepcorredges.dual, pepcorredges.dual.neg, dualmodgenes.vneg, gzallt.cccn.edges, gzallt.cccn.edges.plus, dualpack.vneg, dualpubi.vneg, dualackubi.vneg, gzallt.key, gzallt.gene.key, gzalltgene.data, gzalltgene.ave.ratios, gzgene.cfn.netatts, gzcccn.netatts, gzalltgene.cfn, gzalltgene.cfn.g, gz.cfn, gzallt.all.cf, gzalltgene.all.cf, gz.cf, gzallt.network, file=paste(comp_path, "/Dropbox/_Work/R_/_LINCS/_KarenGuolin/", "TenCell.RData", sep=""))
##############################################################
#_____
which(eu.sp.sed.gzallt.sizes==max(eu.sp.sed.gzallt.sizes)) # 3; 467; 244
head(eu.sp.sed.gzallt.data[[244]])
look <- graph.clust6d.l (clust.data.from.vec(eu.sp.sed.gzallt[[244]], tbl=gzdata.allt))
look <- graph.clust6d.l (clust.data.from.vec(eu.sp.sed.gzallt[[244]], tbl=gzdata.allt))
efractNA.gzall <- sapply(eu.sp.sed.gzallt.data, fractNA)
hist(unlist(efractNA.gzall), breaks=89, col="green3")
fractNA(gzdata.all) # 81% vs 
fractNA(gztencelldata.trimmed)  # .7079
fractNA(gzdata.allt) # 0.77679
# Interesting, 
# Almost all clusters are less than this. Max = 0.7396825
maxone <- which(max(unlist(efractNA.gzall)))
hist(efractNA.gzall, breaks=80, col="maroon")
plot(eu.sp.sed.gzallt.sizes, efractNA.gzall)
which(efractNA.gzall>0.708) # 84
graph.clust6d.l (eu.sp.sed.gzallt.data[which(efractNA.gzall>0.708)]) # H2228 mostly (okay)
which(efractNA.gzall>0.6) 
# 9  84  86  88  90  92 259 261 265 269 295 315 320 339 359 363 420 444 445 448 509 637 710 778
#  All these have at least two columns in common
#contrast:
look <- graph.clust6d.l (eu.sp.sed.gzallt.data[[778]]) # "328.274.205"
look2 <-  graph.clust6d.l(essgzallt.data["328.274.205"])
graph.clust6d.l (eu.sp.sed.gzallt.data[which(efractNA.gzall>0.5)[5]])
which(efractNA.gzall<0.25)
graph.clust6d.l (eu.sp.sed.gzallt.data[which(efractNA.gzall<0.25)[1]])
# The tencell data suggests we should cull clusters with sparse data
# But with gzallt it looks  better
#

fractNA(gzdata.all) # 81% vs 
fractNA(gztencelldata.trimmed)  # .7079
fractNA(gzdata.allt) # 0.77679
# Make more convenient names
# essgzallt <- eu.sp.sed.gzallt # trimmmed above
# Check names, which were lost in proceudre to make the data list
head(eu.sp.sed.gzallt)
lapply(eu.sp.sed.gzallt.data[1:5], rownames) # looks okay, check all
testrownames <- lapply(eu.sp.sed.gzallt.data, rownames)
identical (testrownames, eu.sp.sed.gzallt) # FALSE because names are differet
testcontents <- eu.sp.sed.gzallt
names(testcontents) <- NULL
identical(testrownames, testcontents) # Still FALSE
setdiff(testrownames, testcontents) # Several clusters with NAs
# _____
# neg cor edges and netatts from KGnegcorredges2.R
###################################################################################################
# Compare different results from different data sets/treatments
# "Raw" = without ratios
# Orignial data clusters: eu.sp.sed.gz (306); eu.sp.sed.gzraw (267)
# Ten cellline data clusters: eu.sp.sed.gztencell (612)
# Both data sets clusters: eu.sp.sed.gzallt (839) renamed essgzallt (trimmed to 818)
ratiocommon.1 <- list.common(eu.sp.sed.gz, eu.sp.sed.gzraw)  
# 236 clusters 
hist(sapply(ratiocommon.1, length), breaks=40, col="orange")
tencellcommon <- list.common(eu.sp.sed.gz, eu.sp.sed.gztencell) # only 8
onevstwo <- list.common(eu.sp.sed.gz, essgzallt) # 135
hist(sapply(onevstwo, length), breaks=40, col="chartreuse")
length(intersect(rownames(gzdata), rownames(tencelldata))) 
# 1678 out of 12059 total (9215 trimmed)
# Write these lists as text files using sink(), which opens a connection, and must be closed after the print() 
sink(file=paste("_LINCS/_KarenGuolin/", "ratiocommon.txt", sep=""))
print(ratiocommon.1)
sink()
####################################################################################################
# CCCNs
####################################################################################################
#remove self loops (already done): diag(gzallt.cor) <- NA
#
####>>>> Do first with all data including ratios
# Adjacency matrix:
require(plyr)
gzallt.adj <- rbind.fill.matrix(llply(essgzallt, make.adj.mat))
rownames(gzallt.adj) <- colnames(gzallt.adj)
# dim 7398 7398
gzallt.adj.o <- gzallt.adj[order(rownames(gzallt.adj)), order(colnames(gzallt.adj))]
gzallt.cccn.1 <- gzallt.cor[rownames(gzallt.cor) %in% rownames(gzallt.adj.o), colnames(gzallt.cor) %in% colnames(gzallt.adj.o)]
setdiff(rownames(gzallt.adj), rownames(gzallt.cccn.1)) # 0
length(intersect(rownames(gzallt.adj), rownames(gzallt.cor)))  # 7398 all good
gzallt.cccn <- gzallt.cor[intersect(rownames(gzallt.adj.o), rownames(gzallt.cor)), intersect(colnames(gzallt.adj.o), colnames(gzallt.cor))]
identical(rownames(gzallt.cccn), colnames(gzallt.cccn))
# [1] TRUE
dim(gzallt.cccn)
# 7398 7398
gzallt.NA <- which(is.na(gzallt.adj.o), arr.ind = TRUE)
gzallt.cccn <- replace (gzallt.cccn, gzallt.NA, NA)
# remove self loops
if(any(!is.na(diag(gzallt.cccn)))) {diag(gzallt.cccn) <- NA}
## ****
# Option: Limit to a particlar correlation value like this.
#gzallt.cccn.halflim <- replace (gzallt.cccn, abs(gzallt.cccn)<0.5, NA)
# Make igraph object
gzallt.cccn0 <- gzallt.cccn
  gzallt.cccn0[is.na(gzallt.cccn0)] <- 0
gzallt.cccn.g <- graph.adjacency(as.matrix(gzallt.cccn0), mode="lower", diag=FALSE, weighted="Weight")
#  
hist(edge_attr(gzallt.cccn.g)[[1]], breaks=1000, col="deepskyblue", ylim=c(0,1200))
gzallt.cccn.edges <- 
#-----------------------------------
# 
# Construct CCCN for proteins (genes) based on correlations of mod sites: 
#	
any(duplicated(rownames(gzallt.cccn)))		#  FALSE
any(duplicated(colnames(gzallt.cccn)))		#  FALSE
# Make Gene CCCN
# Note: Did this on server for large data set, much faster at ddply step
gzallt.gene.cccn <- data.frame(gzallt.cccn, row.names = rownames(gzallt.cccn), check.rows=TRUE, check.names=FALSE, fix.empty.names = FALSE)
identical(rownames(gzallt.gene.cccn), colnames(gzallt.gene.cccn)) # TRUE 
# NOTE: if not this can help if check.names=TRUE:
# names(gzallt.gene.cccn) <- sapply(names(gzallt.gene.cccn), function (x) gsub("\\.", " ", x))
# This may be a problem only with the first draft in which not all gene names are correctly identified. A workaround:
# xnames <- names(gzallt.gene.cccn)[startsWith(names(gzallt.gene.cccn),"X")]
# badnames <- xnames[grep("X[[:digit:]]", xnames)]
# goodnames <- substring(badnames, 2, )
# names(gzallt.gene.cccn)[which(names(gzallt.gene.cccn) %in% badnames)] <-  goodnames
# gzallt.gene.cccn <- gzallt.gene.cccn[sort(rownames(gzallt.gene.cccn)), sort(names(gzallt.gene.cccn))]
# identical(rownames(gzallt.gene.cccn), colnames(gzallt.gene.cccn)) # should not be FALSE
# Note in first draft there are many cases of ambiguous gene names separated by semicolons.
gzallt.gene.cccn$Gene.Name <- sapply(rownames(gzallt.gene.cccn), function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])
length(unique(gzallt.gene.cccn$Gene.Name))   # 2913 
any(is.na(gzallt.gene.cccn$Gene.Name)) # F
# gzallt.gene.cccn$Gene.Name[grep("NA", gzallt.gene.cccn$Gene.Name, fixed=TRUE)]
# All genes with NA in the name, e.g. "CTNNA1"
# Use only upper triangle so correlations are not duplicated during the next step
gzallt.gene.cccn[lower.tri(gzallt.gene.cccn)] <- NA	
# ________________________________
# Sum correlations in one dimension, then the other dimension
gzallt.gene.cccn2 <- ddply(gzallt.gene.cccn, .(Gene.Name), numcolwise(function(x) sum(x, na.rm=T)), .progress = "tk")
dim(gzallt.gene.cccn2)	#	2913 7399
any(is.na(gzallt.gene.cccn2$Gene.Name)) # F
rownames(gzallt.gene.cccn2) <- gzallt.gene.cccn2$Gene.Name
gzallt.gene.cccn2 <- gzallt.gene.cccn2[, 2:ncol(gzallt.gene.cccn2)]
gzallt.gene.cccn2 <- data.frame(t(gzallt.gene.cccn2))
gzallt.gene.cccn2$Gene <- sapply(rownames(gzallt.gene.cccn2), function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])
any(is.na(gzallt.gene.cccn2$Gene)) # F
# other dimension
gzallt.gene.cccn3 <- ddply(gzallt.gene.cccn2, .(Gene), numcolwise(function(x) sum(x, na.rm=T)), .progress = "tk")
dim(gzallt.gene.cccn3)  #  2913 2914
any(is.na(gzallt.gene.cccn3$Gene)) # F
# Check row and column names alignment 
print(gzallt.gene.cccn3[1:20, 1:8])
tail(data.frame(gzallt.gene.cccn3$Gene, names(gzallt.gene.cccn3[2:ncol(gzallt.gene.cccn3)])))
# R likes to put dots in column names, which is a problem for ambiguous gene names and gene names with hyphens 
# E.g., Note punctuation missing for "HLA"  "NKX2" "NME1 in column names
gzallt.gene.cccn3$Gene[grep("HLA", gzallt.gene.cccn3$Gene)]
names(gzallt.gene.cccn3)[grep("HLA", names(gzallt.gene.cccn3))]
# ,,,
# Okay this is a pain, so just work around the problem (once satisfied that the gene names actually match). 
names(gzallt.gene.cccn3)[2:ncol(gzallt.gene.cccn3)] <- gzallt.gene.cccn3$Gene
rownames(gzallt.gene.cccn3) <- gzallt.gene.cccn3$Gene
identical (rownames(gzallt.gene.cccn3), names(gzallt.gene.cccn3[,2:ncol(gzallt.gene.cccn3)])) # TRUE
# ________________________________
# Make a adjacency matrix with NAs to index edges only for genes that co-cluster
gzallt.gene.cccn0 <- gzallt.gene.cccn3[,2:ncol(gzallt.gene.cccn3)]
zero.to.NA <- function(df) {
  zer0 <- which(df==0, arr.ind = TRUE)
  cfNA <- as.matrix(df)
  cfNA[zer0] <- NA
  cfNA <- data.frame(cfNA)
  return(cfNA)
}
gzallt.gene.cccn.na <- zero.to.NA(gzallt.gene.cccn0)
# Check
gzallt.gene.cccn.na[1:20, 1:8]
hist(unlist(gzallt.gene.cccn.na), breaks=1000, col="red")
hist(unlist(gzallt.gene.cccn0), breaks=1000, col="blue", ylim=c(0, 1200))
# Make igraph object
gzallt.gene.cccn.g <- graph.adjacency(as.matrix(gzallt.gene.cccn0), mode="lower", diag=FALSE, weighted="Weight")
#  
hist(edge_attr(gzallt.gene.cccn.g)[[1]], breaks=1000, col="deepskyblue", ylim=c(0,1200))
# to get the adjacency matrix back, but NOTE, edges are 0 or 1.
# gzallt.gene.cccn.mat <- as_adjacency_matrix(gzallt.gene.cccn.g, type="both")
# ________________________________
# Here is how to make an edge list file
gzalltgenecccn.edges <- data.frame(as_edgelist(gzallt.gene.cccn.g))
names(gzalltgenecccn.edges) <- c("Gene.1", "Gene.2")
gzalltgenecccn.edges$Weight <- edge_attr(gzallt.gene.cccn.g)[[1]]
gzalltgenecccn.edges$interaction <- "correlation" 
gzalltgenecccn.edges $interaction[gzalltgenecccn.edges$Weight<=-0.5] <- "negative correlation"
gzalltgenecccn.edges $interaction[gzalltgenecccn.edges$Weight>=0.5] <- "positive correlation"
# dim 60282     4
# in above: save(gzalltgenecccn.edges, gzallt.gene.cccn.g, gzallt.gene.cccn, gzallt.cccn.g, gzallt.cccn, )
#######################################################################
# Make CFN with PPI edges
# Test whether combined.all.ppi from KGFunDataObjects.RData covers the genes
gzallt.genes.1 <- unique(c(gzalltgenecccn.edges$Gene.1, gzalltgenecccn.edges$Gene.2))
combined.ppi.genes <- unique(c(combined.all.ppi$Gene.1, combined.all.ppi$Gene.2))
leftover <- setdiff(gzallt.genes.1, combined.ppi.genes)  
# 20 genes not in combined; 4 HLA isoforms, most not in "function.key"

#######################################################################
# Write files for perusal 
write.table(blr.raw, file= "Bin_lung_rawratios.txt", sep="\t", eol = "\n", quote = FALSE, row.names=FALSE, col.names=TRUE)

#---------+++===>>> Plot betweenness to compare different cfns.