# eu.sp.sed.gzallt is the cluster list
# eu.sp.sed.gzallt.data is the cluster list with data. Prune to look at ratios. 
test <- eu.sp.sed.gzallt.data[[2]]
test.mat <- as.matrix(test[, names(test)[grep("Ratio", names(test))]])
# Examine clusters' ratio data: tencellratios.lim.log2.trimmed plus previous data
gzdata.allt.ratios <- gzdata.allt[,names(gzdata.allt)[grep("atio",names(gzdata.allt))]]
esp.gz.ratios <- lapply (eu.sp.sed.gzallt, function (x) gzdata.allt.ratios[rownames(gzdata.allt.ratios) %in% x,])
esp.gz.ratios <- lapply(esp.gz.ratios, function (df) Filter(function(x) !all(is.na(x)), df))
whichclusters <- lapply(esp.gz.ratios, function (x) any(x>=log2(2.25)))
drugchangedratios <- esp.gz.ratios[whichclusters==TRUE]
# 707/839
drugchangedratios <- lapply(drugchangedratios, function (df) Filter(function(x) !all(is.na(x)), df))
# Now filter rows
test=whichdrugs$`402.28.82`
test=drugchangedratios["402.28.82"]
look <- graph.clust6d.l(test)
apply(test, 1, function (x) x>=log2(2.25))
whichptms <- t(apply(test, 1, function (x) abs(x)>=log2(2.25)))
getptms <- which(whichptms, arr.ind = TRUE)
peps <- unique(rownames(getptms))
drugratios <- names(test)[unique(getptms[,2])]
test[peps,drugratios]
look <- graph.clust6d.l(test[peps,drugratios])
# Make this into a function to lapply to list
pick.drug.affected <- function (df, changelimit=log2(2.25)) {
  df <- Filter(function(x) !all(is.na(x)), df)
  if (ncol(df)==0 | all(is.na(df))) {return(NA)} else {
  whichptms <- t(apply(df, 1, function (x) abs(x)>=changelimit))
  getptms <- which(whichptms, arr.ind = TRUE)
  peps <- unique(rownames(getptms))
  if(length(peps)==0) {peps <- names(data.frame(whichptms))}
  drugratios <- names(df)[unique(!is.na(getptms[,2]))]
  return(df[peps,drugratios]) }
}
#drugaffected <- lapply(whichdrugs, pick.drug.affected)
drugaffected <- lapply(esp.gz.ratios, pick.drug.affected)
goodbad <- lapply(drugaffected, function (x) !all(is.na(x)))
drugaffected <- drugaffected[goodbad==TRUE]
# save this to a text file

cat(capture.output(print(drugaffected), file="/Users/_mark_/Dropbox/_Work/R_/_LINCS/_KarenGuolin/drugaffectedPTMs.txt"))
# Which genes are in these clusters?
drugaffectedgenes <- lapply(drugaffected, function (x) get.gene.names.from.peps(rownames(x)))
drugsthataffectgenes <- lapply(drugaffected, names)
cat(capture.output(print(drugaffectedgenes), file="/Users/_mark_/Dropbox/_Work/R_/_LINCS/_KarenGuolin/drugaffectedgenes.txt"))
cat(capture.output(print(drugsthataffectgenes), file="/Users/_mark_/Dropbox/_Work/R_/_LINCS/_KarenGuolin/drugsthataffectPTMs.txt"))
#-----------------------
# What can we say about PTMs and drugs that affect different PTMs? 
# >>>>>>>>>>>
afatcols <- names(gzdata.allt.ratios)[grep("Afat", names(gzdata.allt.ratios))]
crizcols <- names(gzdata.allt.ratios)[c(1:3, grep("Criz", names(gzdata.allt.ratios)))]
erlotcols <- names(gzdata.allt.ratios)[c(7:9, grep("Erlot", names(gzdata.allt.ratios)))]
PRcols <- names(gzdata.allt.ratios)[grep("SEPTM.PR", names(gzdata.allt.ratios))]
dasatcols <- names(gzdata.allt.ratios)[grep("Dasat", names(gzdata.allt.ratios))]
# check
names(gzdata.allt.ratios) %w/o% c(afatcols, crizcols, dasatcols, erlotcols,  PRcols) # 0
#
# What is the intersection of drug-affected PTMs in clusters?
test.d <- (drugaffected$`7.29.7`)
dim(test.d[, names(test.d) %in% c(erlotcols, crizcols)]) # 10/13
test.d.ec <- test.d[, names(test.d) %in% c(erlotcols, crizcols)]
test.d.eandc <- test.d[, names(test.d) %in% erlotcols & names(test.d) %in% crizcols]
test.d.e <- test.d[, names(test.d) %in% erlotcols]
test.d.c <- test.d[, names(test.d) %in% crizcols]
e.ptms <- rownames(test.d.e)[which(apply(test.d.e, 1, function (x) !all(is.na(x))))]
d.ptms <- rownames(test.d.c)[which(apply(test.d.e, 1, function (x) !all(is.na(x))))]
intersect(e.ptms, d.ptms) # all of them
outersect(e.ptms, d.ptms) # 0

# Find them as above in all data
whichptms <- t(apply(gzdata.allt.ratios, 1, function (x) abs(x)>=log2(2.25)))
getptmdata <- which(whichptms, arr.ind = TRUE)
peps <- unique(rownames(getptmdata))
drugs <- names(gzdata.allt.ratios)[unique(getptmdata[,2])]
getptmdata <- data.frame(getptmdata)
getptmdata$drugs <- names(gzdata.allt.ratios)[getptmdata$col]
getptmdata$afatinib <- getptmdata$drugs %in% afatcols
getptmdata$crizotinib <- getptmdata$drugs %in% crizcols
getptmdata$dasatinib <- getptmdata$drugs %in% dasatcols
getptmdata$erlotinib <- getptmdata$drugs %in% erlotcols
getptmdata$PR171 <- getptmdata$drugs %in% PRcols
#
head(gzdata.allt.ratios[getptmdata[1,1], getptmdata[1,2]])
change <- vector()
for (i in 1:dim(getptmdata)[1]){
  change[i] <- gzdata.allt.ratios[getptmdata[i,1], getptmdata[i,2]]
}
getptmdata$change <- change
getptmdata$phosphorylation <- grepl(".p.", rownames(getptmdata), fixed=TRUE)
getptmdata$acetylation <- grepl(".ack.", rownames(getptmdata), fixed=TRUE)
getptmdata$ubiquitination <- grepl(".ubi.", rownames(getptmdata), fixed=TRUE)
# preserve rownames 
rownames(getptmdata)[grep("Y440", rownames(getptmdata))] # NOTE! ambigous names!
# getptmdata$Peptide.Name <- gsub("\\.", " ", rownames(getptmdata), fixed=F) 
# leaves trailing numbers - 
# have to get rid of terminal numbers (from duplicates) first
# Note: this will also remove ambigous names
Peptide.Name.list <- strsplit(rownames(getptmdata), "\\.")
# Check ambigous names
Peptide.Name.list[grep("Y440", Peptide.Name.list)]
# Peptide.Name <- lapply(Peptide.Name.list, function (x) paste(x[1:3], collapse=" " ))
Peptide.Name[grep("Y440", Peptide.Name)]
# NOTE: paste(x[1:3] removes ambigous names
Peptide.Name <- lapply(Peptide.Name.list, function (x) paste(x[1:length(x)], collapse=" " ))
Peptide.Name[grep("Y440", Peptide.Name)]
Peptide.Name.list[grep("Y440", Peptide.Name.list)]
# Note that the ";" is now gone
# The gsub/str_replace methods are not working, try searching
gzdata.allt.ratios[grep("FYN p Y440", rownames(gzdata.allt.ratios)), 1:2]
rownames(gzdata.allt.ratios[grep("FYN p Y440", rownames(gzdata.allt.ratios)), ])
return.ptm <- function(element) {
  el.length <- length(element)
  pepname <- paste(element[1:3], collapse=" " )
  if (el.length>3) {
    pepname <- rownames(gzdata.allt.ratios[grep(pepname, rownames(gzdata.allt.ratios)), ])
  }
  return(pepname)
}
Peptide.Name <- lapply(Peptide.Name.list, return.ptm)
Peptide.Name[grep("Y440", Peptide.Name)]  # Works!!


#>>>>>
Peptide.Name[which(grepl("NA ", Peptide.Name))]
whichNA <- which(grepl("NA ", Peptide.Name))
getptmdata[whichNA,]
# Why are there NAs? 
whichrowNA <- which(grepl("NA ", rownames(gzdata.allt)))
rownames(gzdata.allt)[whichrowNA]            
# 10 with NA in names.                     
badnames <- c("NA ubi K162; NA ubi K132", "NA ubi K359", "NA ubi K252; NA ubi K338; TUBB4A ubi K252; TUBB ubi K252; TUBB4B ubi K252; TUBB3 ubi K252; TUBB2A ubi K252; TUBB8 ubi K252; TUBB6 ubi K252; TUBB2B ubi K252")  
morebadnames <- rownames(gzdata.allt)[whichrowNA[c(41:50)]]
  # NOTE: Remvoe these XXXX                  
getptmdata$Peptide.Name <- Peptide.Name

getptmdata.df <- getptmdata[-(getptmdata$Peptide.Name %in% c(badnames, morebadnames)),]
# 
# Check
all.tki.ptms.df <- getptmdata.df[which(apply(getptmdata.df[,4:7], 1, all)),]
# !!! No single PTM is in all columns
ptmdata <- dlply(getptmdata, .(afatinib, crizotinib, dasatinib, erlotinib, PR171, phosphorylation, acetylation, ubiquitination))
head(ptmdata[[1]]["change"])
sortnames <- names(getptmdata)[c(4:8, 10:12)]
combos <- apply(getptmdata[,sortnames], 1, function(x) length(which(x))) # all 2 okay
datanames <- c(noquote(strsplit(names(ptmdata), "\\.")))
sortnames[as.logical(datanames[[1]])]
newnames <- lapply(datanames, function (x) sortnames[as.logical(x)])
names(ptmdata) <- newnames
ptmchanges <- sapply(ptmdata, "[", 'change')
names(ptmchanges) <- lapply(newnames, function(x) paste(noquote(x), collapse=", "))

## What is the overlap of PTMs affected by drugs?
# combinations: 
phosnames <- c(3, 6, 9, 12, 15)
acknames <- c(2, 5, 8, 11, 14)
ubnames <- c(1, 4, 7, 10, 13)
ptms.by.drug <- sapply(ptmdata, "[", 'Peptide.Name')
names(ptms.by.drug) <- names(ptmchanges)
# This functions computes intersect of all list elements between two lists
commonpeps.p1 <- list.common(ptms.by.drug[phosnames], ptms.by.drug[phosnames], keeplength = 1)
sapply(commonpeps.p1, length)
sapply(ptms.by.drug[phosnames], length)
sapply(ptms.by.drug, length) 
# **** Sum for paper
phossites <- sapply(ptms.by.drug[phosnames], length) 
sum(phossites[2:5])
# 4987
acksites <- sapply(ptms.by.drug[acknames], length) 
sum(acksites[2:5])
# 3249
ubsites <- sapply(ptms.by.drug[ubnames], length) 
sum(ubsites[2:5])
# 4452
# Remove comparisons to same list
dupes <- names(commonpeps.p1)[c(1, 7, 13, 19, 25)]
commonpeps.p <- commonpeps.p1[names(commonpeps.p1) %w/o% dupes] 
commonpeps.pp <- lapply(commonpeps.p, unlist)
# Do this again to find interesction of intersections
commonpeps.pp12 <- list.common(commonpeps.pp, commonpeps.pp, keeplength = 1)
hist(sapply(commonpeps.pp12, length), col="magenta", breaks=110)
# now a list of 400! Unweildy
lengths.pp12 <- sapply(commonpeps.pp12, length)
# commonpepspp12.pruned <- commonpeps.pp12[which(lengths.pp12>0)] # no diff
range(lengths.pp12)    # 18 677
# Focus on two drugs
crizerlot.p1 <- unlist(intersect(ptms.by.drug[["crizotinib, phosphorylation"]], ptms.by.drug[["erlotinib, phosphorylation"]]))
crizerlot.a1 <- unlist(intersect(ptms.by.drug[["crizotinib, acetylation"]], ptms.by.drug[["erlotinib, acetylation"]]))
crizerlot.u1 <- unlist(intersect(ptms.by.drug[["crizotinib, ubiquitination"]], ptms.by.drug[["erlotinib, ubiquitination"]]))
# Add a third
dasaterlot.p1 <- unlist(intersect(ptms.by.drug[["dasatinib, phosphorylation"]], ptms.by.drug[["erlotinib, phosphorylation"]]))
dasaterlot.a1 <- unlist(intersect(ptms.by.drug[["dasatinib, acetylation"]], ptms.by.drug[["erlotinib, acetylation"]]))
dasaterlot.u1 <- unlist(intersect(ptms.by.drug[["dasatinib, ubiquitination"]], ptms.by.drug[["erlotinib, ubiquitination"]]))
#
dasatcrizot.p1 <- unlist(intersect(ptms.by.drug[["dasatinib, phosphorylation"]], ptms.by.drug[["crizotinib, phosphorylation"]]))
dasatcrizot.a1 <- unlist(intersect(ptms.by.drug[["dasatinib, acetylation"]], ptms.by.drug[["crizotinib, acetylation"]]))
dasatcrizot.u1 <- unlist(intersect(ptms.by.drug[["dasatinib, ubiquitination"]], ptms.by.drug[["crizotinib, ubiquitination"]]))
# and a fourth
afatcrizot.p1 <- unlist(intersect(ptms.by.drug[["afatinib, phosphorylation"]], ptms.by.drug[["crizotinib, phosphorylation"]]))
afatcrizot.a1 <- unlist(intersect(ptms.by.drug[["afatinib, acetylation"]], ptms.by.drug[["crizotinib, acetylation"]]))
afatcrizot.u1 <- unlist(intersect(ptms.by.drug[["afatinib, ubiquitination"]], ptms.by.drug[["crizotinib, ubiquitination"]]))
#
afatdasat.p1 <- unlist(intersect(ptms.by.drug[["afatinib, phosphorylation"]], ptms.by.drug[["dasatinib, phosphorylation"]]))
afatdasat.a1 <- unlist(intersect(ptms.by.drug[["afatinib, acetylation"]], ptms.by.drug[["dasatinib, acetylation"]]))
afatdasat.u1 <- unlist(intersect(ptms.by.drug[["afatinib, ubiquitination"]], ptms.by.drug[["dasatinib, ubiquitination"]]))
#
afaterlot.p1 <- unlist(intersect(ptms.by.drug[["afatinib, phosphorylation"]], ptms.by.drug[["erlotinib, phosphorylation"]]))
afaterlot.a1 <- unlist(intersect(ptms.by.drug[["afatinib, acetylation"]], ptms.by.drug[["erlotinib, acetylation"]]))
afaterlot.u1 <- unlist(intersect(ptms.by.drug[["afatinib, ubiquitination"]], ptms.by.drug[["erlotinib, ubiquitination"]]))
#___________________
# 1. criz
# 2. erlot
# 3. dasat
# 4. afat
cp <- unlist(ptms.by.drug[["crizotinib, phosphorylation"]])
ep <- unlist(ptms.by.drug[["erlotinib, phosphorylation"]])
dp <- unlist(ptms.by.drug[["dasatinib, phosphorylation"]])
ap <- unlist(ptms.by.drug[["afatinib, phosphorylation"]])
#
library("VennDiagram")
n12 <- unlist(crizerlot.p1)
n13 <- unlist(dasatcrizot.p1)
n14 <- unlist(afatcrizot.p1)
n23 <- unlist(dasaterlot.p1)
n24 <- unlist(afaterlot.p1 )
n34 <- unlist(afatdasat.p1 )

n123 = intersect(n12, dp)
n124 = intersect(n12, ap)
n134 = intersect(n13, ap)
n234 = intersect(n23, ap)

n1234 = intersect(n123, n234)
# Check
all.drug.phos1 <- gzdata.allt.ratios[rownames(gzdata.allt.ratios) %in% n1234,] # only 53? # Fixed: now 60
all.drug.phos <- gzdata.all[gzdata.all$Row.names %in% n1234,] # same
out.1 <- outersect(n1234, all.drug.phos$Row.names)
# now zero
# gzdata.all[gzdata.all$Row.names %in% out.1, 1:4]
cccn.nodes <- unique(c(gzallt.cccnplus$source, gzallt.cccnplus$target)) # 7715
# cccn.nodes[cccn.nodes %in% out.1] # 0 

write.table(all.drug.phos, file="AllDrugAffectedPhos.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")


par(mar=c(4.2, 6.1, 4.1, 6.1), mai=c(2.02, 2.2, 0.82, 0.42), oma=c(0,0,0,0))
dev.new()

vp = draw.quad.venn(length(cp), length(ep), length(dp),  length(ap),  length(n12), length(n13), length(n14), length(n23), length(n24), length(n34), length(n123), length(n124), length(n134), length(n234), length(n1234), category=c("crizotinib","erlotinib", "dasatinib", "afatinib"), 	fill = c("red", "green", "cyan", "lightblue"), fontface = "bold", fontfamily="", cat.fontfamily = rep("", 4), cat.cex = 1.8, cex=2, margin=0.08)

# vp = draw.quad.venn(length(cp), length(ep), length(dp),  length(ap),  length(n12), length(n13), length(n14), length(n23), length(n24), length(n34), length(n123), length(n124), length(n134), length(n234), length(n1234), category=c("crizotinib","dasatinib","erlotinib", "afatinib"), 	fill = c("red", "green", "cyan", "lightblue"), fontface = "bold", fontfamily="", cat.fontfamily = rep("", 4), cat.cex = 1.8, cex=2, margin=0.08)

#___________________
# Acetylation
# 1. criz
# 2. erlot
# 3. dasat
# 4. afat
ca <- ptms.by.drug[["crizotinib, acetylation"]]
ea <- ptms.by.drug[["erlotinib, acetylation"]]
da <- ptms.by.drug[["dasatinib, acetylation"]]
aa <- ptms.by.drug[["afatinib, acetylation"]]
#
n12 <- unlist(crizerlot.a1)
n13 <- unlist(dasatcrizot.a1)
n14 <- unlist(afatcrizot.a1)
n23 <- unlist(dasaterlot.a1)
n24 <- unlist(afaterlot.a1 )
n34 <- unlist(afatdasat.a1 )

n123 = intersect(n12, unlist(dasatcrizot.a1))
n124 = intersect(n12, unlist(afatcrizot.a1))
n134 = intersect(n13, unlist(afatcrizot.a1))
n234 = intersect(n23, unlist(afatcrizot.a1))

n1234 = intersect(n123, n234)
par(mar=c(4.2, 6.1, 4.1, 6.1), mai=c(2.02, 2.2, 0.82, 0.42), oma=c(0,0,0,0))
dev.new()

vp = draw.quad.venn(length(ca), length(ea), length(da),  length(aa),  length(n12), length(n13), length(n14), length(n23), length(n24), length(n34), length(n123), length(n124), length(n134), length(n234), length(n1234), category=c("crizotinib","erlotinib", "dasatinib", "afatinib"), 	fill = c("red", "green", "cyan", "lightblue"), fontface = "bold", fontfamily="", cat.fontfamily = rep("", 4), cat.cex = 1.8, cex=2, margin=0.08)

all.drug.ac <- gzdata.all[gzdata.all$Row.names %in% n1234,]
write.table(all.drug.ac, file="AllDrugAffectedAcetyl.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

#___________________
#___________________
# Ubiqutination
# 1. criz
# 2. erlot
# 3. dasat
# 4. afat
cu <- ptms.by.drug[["crizotinib, ubiquitination"]]
eu <- ptms.by.drug[["erlotinib, ubiquitination"]]
du <- ptms.by.drug[["dasatinib, ubiquitination"]]
au <- ptms.by.drug[["afatinib, ubiquitination"]]
#
n12 <- unlist(crizerlot.u1)
n13 <- unlist(dasatcrizot.u1)
n14 <- unlist(afatcrizot.u1)
n23 <- unlist(dasaterlot.u1)
n24 <- unlist(afaterlot.u1 )
n34 <- unlist(afatdasat.u1 )

n123 = intersect(n12, unlist(dasatcrizot.u1))
n124 = intersect(n12, unlist(afatcrizot.u1))
n134 = intersect(n13, unlist(afatcrizot.u1))
n234 = intersect(n23, unlist(afatcrizot.u1))

n1234 = intersect(n123, n234)
par(mar=c(4.2, 6.1, 4.1, 6.1), mai=c(2.02, 2.2, 0.82, 0.42), oma=c(0,0,0,0))
dev.new()

vp = draw.quad.venn(length(cu), length(eu), length(du),  length(au),  length(n12), length(n13), length(n14), length(n23), length(n24), length(n34), length(n123), length(n124), length(n134), length(n234), length(n1234), category=c("crizotinib","erlotinib", "dasatinib","afatinib"), 	fill = c("red", "green", "cyan", "lightblue"), fontface = "bold", fontfamily="", cat.fontfamily = rep("", 4), cat.cex = 1.8, cex=2, margin=0.08)
#________________

all.drug.ub <- gzdata.all[gzdata.all$Row.names %in% n1234,]
write.table(all.drug.ub, file="AllDrugAffectedUbiquitination.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
# For paper table
all.drug.ptms <- rbind(all.drug.phos, all.drug.ac, all.drug.ub)
write.table(all.drug.ptms , file="AllDrugAffectedPTMs.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

# For supplemental information
write.table(gzdata.allt , file="AlldataPTMs.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
# Graph this?

all.drug.genes <- extract.genes.from.peplist(rownames(all.drug.ptms))
# open CCCN_CFN_Plus_Subs1.cys
selectNodes(all.drug.ptms$Row.names, by="id", preserve=FALSE)
selectFirstNeighbors()
createSubnetwork(nodes = getSelectedNodes(), nodes.by.col = "id", subnetwork.name="PTMs affected by all drugs")
# Interesting selection of more than 20 cliques!

drugcols <- names(getptmdata)[4:8]
twodrugs <- apply(getptmdata[,drugcols], 1, function(x) length(which(x)))
  
  getptmdata[apply(getptmdata[,TFcols], 1, function(x) length(which(c(x))))>1,]


# now make some graphs
hist(getptmdata[getptmdata$phosphorylation, "change"])
hist(getptmdata[(getptmdata$phosphorylation & getptmdata$crizotinib), "change"], breaks=100)
# Grouped data for boxplot

dev.new()
par(mar=c(6.2, 4.1, 4.1, 2.1), mai=c(2.02, 1.2, 0.82, 0.42))
#boxplot(ptmchanges)
#stripchart(ptmchanges)
#dev.new()
plot(density(ptmchanges[[1]]), lwd=4, xlab=names(ptmchanges[1]), cex.lab=1.6, cex.axis=1.6, main="")
abline(v=-0.5, lty=2)
abline(v=0.5, lty=2)
polygon(density(ptmchanges[[1]]), col=alpha("blue", 0.5), border="blue")     
# **** repeat length(ptmchanges) = 15
polygon(density(dualpubi$Weight), col=alpha("blue", 0.5), border="black")     
for (i in 1:length(ptmchanges)){
  dev.new()
  par(mar=c(6.2, 4.1, 4.1, 2.1), mai=c(2.02, 1.2, 0.82, 0.42))
  plot(density(ptmchanges[[i]]), lwd=4, xlab=names(ptmchanges[i]), cex.lab=1.6, cex.axis=1.6, main="")
  abline(v=-log2(2.25), lty=2)
  abline(v=log2(2.25), lty=2)
  # manual toggle:
  # polygon(density(ptmchanges[[i]]), col=alpha("blue", 0.5), border="black")     
  # polygon(density(ptmchanges[[i]]), col="yellow", border="black")     
}
polygon(density(ptmchanges[[i]]), col="yellow", border="black")     
polygon(density(ptmchanges[[i]]), col=rainbow(10), border="black")     
polygon(density(ptmchanges[[i]][which(ptmchanges[[i]]<(-log2(2.25)))]), col=alpha("blue", 0.5), border="black", add=T)     
polygon(density(ptmchanges[[i]][which(ptmchanges[[i]]>(log2(2.25)))]), col=alpha("yellow", 0.5), border="black", add=T)     
polygon(c(density(ptmchanges[[i]][which(ptmchanges[[i]]<(-log2(2.25)))]), col=alpha("blue", 0.5), border="black", density(ptmchanges[[i]][which(ptmchanges[[i]]>(log2(2.25)))]), col=alpha("yellow", 0.5), border="black"))    
# Try this
for (i in 1:length(ptmchanges)){
  dev.new()
  par(mar=c(6.2, 4.1, 4.1, 2.1), mai=c(2.02, 1.2, 0.82, 0.42))
  dens <- density(ptmchanges[[i]])
  plot(dens, lwd=4, xlab=names(ptmchanges[i]), cex.lab=1.6, cex.axis=1.6, main="")
  # q2     <- max(ptmchanges[[i]])
  q2 <- max(dens$x) - 0.1
  q65    <- log2(2.25)
  qn08   <- -log2(2.25)
  # qn02   <- min(ptmchanges[[i]])
  qn02 <- min(dens$x) + 0.1

  x1 <- min(which(dens$x >= q2))  
  x2 <- max(which(dens$x <  q65))
  x3 <- min(which(dens$x >= qn08))  
  x4 <- max(which(dens$x <  qn02))

  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="yellow", border="transparent"))
  with(dens, polygon(x=c(x[c(x3,x3:x4,x4)]), y= c(0, y[x3:x4], 0), col=alpha("blue", 0.5), border="transparent"))
  }
#
# Stripchart and overlayed boxplots for drug effects in aggregate
#stripchart(ptmchanges, method="jitter", vertical=TRUE, xlab=NULL, pch=20, col=alpha("grey80", 0.2), jitter=1/3)
stripchart(ptmchanges, method="jitter", vertical=TRUE, xaxt="n", xlab=NULL, ylab="log2 ratio", cex.lab=1.2,  pch=20, col=alpha("blue", 0.1), jitter=1/3)
text(1:15, par("usr")[3] - 0.15, srt=45, adj=c(1, 1),  labels=names(ptmchanges), xpd=T, cex=1.2)
# Boxplot separate positive and negative changes
par(mar=c(6.2, 4.1, 4.1, 2.1), mai=c(2.02, 1.2, 0.82, 0.42), bg="grey", fg="black")

ptmposchanges <- lapply(ptmchanges, function(x) x[which(x>=log2(2.25))])
ptmnegchanges <- lapply(ptmchanges, function(x) x[which(x<=log2(2.25))])
stripchart(ptmchanges, method="jitter", vertical=TRUE, xaxt="n", xlab=NULL, ylab="log2 ratio", cex.lab=1.2,  col=alpha("grey80", 0.2), jitter=1/3)
stripchart(ptmposchanges, method="jitter", vertical=TRUE, xaxt="n", cex.lab=1.2, xlab=NULL,  cex.lab=1.2,  pch=20, col=alpha("yellow", 0.2), jitter=1/3, add=T)
stripchart(ptmnegchanges, method="jitter", vertical=TRUE, xaxt="n", xlab=NULL, cex.lab=1.2,  pch=20, col=alpha("blue", 0.1), jitter=1/3, add=T)
# add boxplots
boxplot(ptmposchanges, add=TRUE,  xaxt="n", col="transparent", fg="red")
boxplot(ptmnegchanges, add=TRUE,  xaxt="n", col="transparent")
text(1:15, par("usr")[3] - 0.15, srt=45, adj=c(1, 1),  labels=names(ptmchanges), xpd=T, cex=1.2)



stripchart(dualtest3.list, method="jitter", vertical=TRUE, pch=20, col=alpha("blue", 0.2), jitter=1/3, add=TRUE)
#stripchart(dualtest4.list, method="jitter", vertical=TRUE, pch=20, col=alpha("blue", 0.2), jitter=1/3, add=TRUE)
text(1:3, par("usr")[3] - 0.15, srt=45, adj=c(1, 1),  labels=c("me me", "me ac", "ac ac"), xpd=T, cex=1.2)
# no.nodes as circle size
symbols(eprf.cfn.tbl$thresholds, eprf.cfn.tbl$density,  circles=eprf.cfn.tbl$no.nodes/40000, inches = FALSE, bg=adjustcolor("green", 0.6), fg = adjustcolor("forestgreen", 0.9), add=TRUE)
title("Cluster Filtered Networks from Ratio Data")


#-----------------------
# Now look at interesction with pathways
keggcommon <- clust.common(drugaffectedgenes, keggpath)
# Error, try it manually
list1=drugaffectedgenes
list2=wikipath
# repeat below using 
list2=keggpath
#list.common <- function (list1, list2, keeplength=3) {
  parse <- lapply(list1, function (y) sapply(list2,  function(x) intersect(x, y)))
  dims <- lapply(parse, function (x) sapply(x, length))
  keep <- which(sapply(dims, sum) > keeplength)
  pare <- parse[keep]
  prune <- lapply(pare, function (y) return (y[which(sapply(y, function (x) which(length(x) > keeplength )) > 0)]))
  newlist <- unlist(prune, recursive=FALSE)
#  return(newlist) }

wikicommon=newlist
cat(capture.output(print(wikicommon), file="/Users/_mark_/Dropbox/_Work/R_/_LINCS/_KarenGuolin/drugaffectedWikiPathGenes.txt"))
keggcommon=newlist
cat(capture.output(print(keggcommon), file="/Users/_mark_/Dropbox/_Work/R_/_LINCS/_KarenGuolin/drugaffectedKEGGPathGenes.txt"))
#-----------------------------
# Which PTMs are affected by drugs and have negative correlations with other PTMs in the same protein? 
gzneg.edges <- rbind(dualpack.vneg[,1:4], dualpubi.vneg[,1:4], dualackubi.vneg[,1:4])
gzneg.ptms <- unique(c(gzneg.edges[,1], gzneg.edges[,2]))
drug.neg.ptms <- lapply(drugaffected, function (y) intersect(rownames(y), gzneg.ptms))
dnp.lengths <- lapply(drug.neg.ptms, length)
drug.neg.ptms  <- drug.neg.ptms [which(dnp.lengths > 0)]
cat(capture.output(print(drug.neg.ptms), file="/Users/_mark_/Dropbox/_Work/R_/_LINCS/_KarenGuolin/drugaffectednegcorrPTMs.txt"))

save(gzdata.allt.ratios, esp.gz.ratios, drugchangedratios, pick.drug.affected, drugaffected, drugaffectedgenes, drugsthataffectgenes, wikicommon, keggcommon, gzneg.edges, gzneg.ptms, drug.neg.ptms, afatcrizot.p1, afatdasat.p1, afaterlot.p1, crizerlot.p1, dasatcrizot.p1, dasaterlot.p1, afatcrizot.a1, afatdasat.a1, afaterlot.a1, crizerlot.a1, dasatcrizot.a1, dasaterlot.a1, afatcrizot.u1, afatdasat.u1, afaterlot.u1, crizerlot.u1, dasatcrizot.u1, dasaterlot.u1, file="/Users/_mark_/Dropbox/_Work/R_/_LINCS/_KarenGuolin/drug_effects.RData")

graph.clust6d.l(drugaffected$`130.57.111`)
# 
graph.clust6d.l(drugaffected$`130.141.111`)
graph.clust6d.l(drugaffected$`55.62.56`)
graph.clust6d.l(drugaffected$`7.29.7`)
graph.clust6d.l(drugaffected$`7.29.7`[drug.neg.ptms$`7.29.7`,] )
graph.clust6d.l(drugaffected$`55.62.56`[drug.neg.ptms$`55.62.56`,] )
# Examine ALL clusters for drug affected PTMs
look <- list()
for (i in 1:length(drugaffected)) {
  if(dim(drugaffected[[i]][1])>2 & dim(drugaffected[[i]][2])>2) {
  graph.clust6d.l(drugaffected[[i]]) }
}
#-----------------------------

# Explore FASN links in Nuclear Receptors and EGFR.cys
  fasnratios <- gzdata.allt.ratios[grep("FASN ", rownames(gzdata.allt.ratios), fixed=TRUE),]
  look <-  graph.clust6d.l(fasnratios)
  lookprune <- look[apply(look, 1, function(x) !all(is.na(x))),]
  lookp <-  graph.clust6d.l(lookprune)
  lookp2 <- look[apply(lookp, 1, function(x) (filled(x))>2),]
  lookp <-  graph.clust6nodnosort.l(lookp2)
  fasnnodes <- getAllNodes()
  fasnedges <- getTableColumns('edge', c("source", "target", "Weight", "interaction"))
  fasnallnegedges <- filter.edges.0(rownames(fasnratios), gzneg.edges)
  lookp3 <-  graph.clust6nodnosort.l(lookp[rownames(lookp) %in% c(fasnallnegedges[,1], fasnallnegedges[,2]),])
  lookp4 <- graph.clust6nodnosort.l(lookp3[apply(lookp3, 2, function(x)  (filled(x))>2),])
  look1 <- graph.clust6d.l(fasnratios[rownames(fasnratios) %in% c(fasnallnegedges[,1], fasnallnegedges[,2]),])
  look2 <- graph.clust6d.l(fasnratios[rownames(fasnratios) %in% fasnnodes,])
  
  
#-----------------------------
# Questions that come to mind with an eye towards experimental approaches that are testable.
  # Are there sets genes that are affected by one drug, and complementary genes in the same pathway that are affected by a different drug, that suggests a two-drug combination therapy? For example, if we treat cells with erlotinib and titrate smaller to larger doses of dasatinib, can we downregulate a paritcular pathway in a synergistic fasion as predicted by the intersection with PTM clusters and pathway genes? 
  # Similarly, for pathways that involve the HSP90 genes, can we predict a synergistic efect of a drug that affects these and titration of smaller to larger doses of geldanamycin (or another HSP90 inhibitor)?
  # Can we predict based on the data synergistic drug combinations that would maximally impact PTMs and activity of signal transduction control points? Similarly, impact FASN and EP300? 
  # The PI3K-Akt Signaling Pathway would be expected to affect cell growth rate. Can we draw pathways to, for example clusters like $`255.57.111 in drugaffectedWikiPathGenes.txt that might help identify synergistic drug combinations to decrease cell growth rate?

  