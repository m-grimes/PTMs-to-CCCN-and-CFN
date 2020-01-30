Data Input and Formatting
================
Mark Grimes
1/30/2020

``` r
#  use install.packages("package_name") to download the package
#  or, for bioconductor packages:
#  To install core packages, type the following in an R command window:
  if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
```

    ## Loading required namespace: BiocManager

``` r
BiocManager::install()
```

    ## Bioconductor version 3.9 (BiocManager 1.30.10), R 3.6.2 (2019-12-12)

``` r
#Install specific packages, e.g., “GenomicFeatures” and “AnnotationDbi”, with
# BiocManager::install(c("Hmisc", "RCy3", "impute", "PMA", "GO.db", "Category", "GOstats"))
# library ("tidyverse") has the following packages embedded 
    #•  ggplot2, for data visualisation.
    #•  dplyr, for data manipulation.
    #•  tidyr, for data tidying.
    #•  readr, for data import.
    #•  purrr, for functional programming.
    #•  tibble, for tibbles, a modern re-imagining of data frames.
    #•  stringr, for strings.
    #•  forcats, for factors.   
library(plyr)
require(gplots)
```

    ## Loading required package: gplots

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::arrange()   masks plyr::arrange()
    ## x purrr::compact()   masks plyr::compact()
    ## x dplyr::count()     masks plyr::count()
    ## x dplyr::failwith()  masks plyr::failwith()
    ## x dplyr::filter()    masks stats::filter()
    ## x dplyr::id()        masks plyr::id()
    ## x dplyr::lag()       masks stats::lag()
    ## x dplyr::mutate()    masks plyr::mutate()
    ## x dplyr::rename()    masks plyr::rename()
    ## x dplyr::summarise() masks plyr::summarise()
    ## x dplyr::summarize() masks plyr::summarize()

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
    "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
    without <- function(x, y) x[!x %in% y] #--  x without y
    nmissing <- function(x) sum(is.na(x))
    filled <- function (x) {length(x) - nmissing(x)}
    fractNA <- function(df) {
         result <- nmissing(df)/(dim(df)[1]*dim(df)[2])
         return(result)
    }
    mean.na <- function(x) mean(x, na.rm=TRUE)
    max.na <- function(x) max(x, na.rm=TRUE)
    min.na <- function(x) min(x, na.rm=TRUE)
     sd.na <- function(x) sd(x, na.rm=TRUE)
  outersect <- function(x,y){sort(c(setdiff(x,y), setdiff(y,x)))}
  
zero.to.NA <- function(df) {
  zer0 <- which(df==0, arr.ind = TRUE)
  cfNA <- as.matrix(df)
  cfNA[zer0] <- NA
  cfNA <- data.frame(cfNA)
  return(cfNA)
}
```

``` r
# requires data table with rownames 
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
        return (clust.data) }
    # end clust.data.from.vec
#+
# A general version of get.gene.name; 
# use with sapply   
# Fixes some common misnamed genes and the problem that Excel turns some gene names into dates.
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
# Make peptide names using this function (revised to Karen's reformatting). 
# Example :
# tencellack.head$Peptide.Name <- mapply(name.peptide, genes=tencellack.head$LeadingGeneSymbols, sites= tencellack.head$Positions, modification="ack", aa=tencellack.head$Amino.acid)
name.peptide <- function (genes, modification="p", sites, aa)   {
  genes.v <- unlist(strsplit(genes, ";", fixed = TRUE))
  sites.v <- unlist(strsplit(sites, ";", fixed = TRUE))
  sites.v <- sapply(sites.v, function (x) paste (aa, x, sep=""))
  Peptide.v <- as.character(noquote(paste(genes.v[1:length(genes.v)], modification, sites.v[1:length(sites.v)], sep=" ")))
  Peptide <- paste(unique(Peptide.v), collapse="; ")
  return(Peptide)
}
```

## Heatmap graphing functions
