rm(list=ls())
library(plyr)
library(dplyr)
library(igraph)

localpath <- "~/Dropbox/Mark/Code/GuolinData"
setwd(localpath)
source("Functions_111218.R")
load(file=paste(localpath,"NewKGData.Rdata", sep="/"))
#load(file=paste(localpath, "LC_TMT_Functions.RData", sep="/"))
#kglog2data -- cleaned data file to work with

#create a directory to save plots
plots_path_dir <- file.path(getwd(), "Plots")
if(!dir.exists(plots_path_dir)){ dir.create(plots_path_dir)}

tsne.res <- calc.tsne(ptm.data=kglog2data, plots_path_dir = plots_path_dir,
          dims = 3, perplexity = 15, theta = 0.25, 
          check_duplicates = FALSE, pca=FALSE, toolong = 3.5)
