Networks: CCCN and CFN Generation
================
Mark Grimes
1/30/2020

## R Markdown

``` r
# New versions of MGRCyFunctions.R linking RCy3 to Cytoscape
# author: Mark Grimes
# This is the best place to report issues:
# https://github.com/cytoscape/RCy3/issues
# Thanks!    - Alex
# Swagger: http://localhost:1234/v1/swaggerUI/swagger-ui/index.html?url=http://localhost:1234/v1/swagger.json#/
# ------------------------------------------------------------------------
# To UpDATE
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("RCy3")
# biocLite("BiocStyle")
# Try the latest version of RCy3 by running:
    # install.packages("BiocManager")
    #   BiocManager::install("RCy3")
library(devtools)
```

    ## Loading required package: usethis

``` r
library(RCy3)
library(plyr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(RColorBrewer)
library(gplots)
```

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
library(igraph)
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library(knitr)
library("BiocStyle")
options(stringsAsFactors=FALSE)
#
# There are several resources available:
browseVignettes('RCy3')
```

    ## starting httpd help server ...

    ##  done

``` r
# https://github.com/cytoscape/RCy3/wiki/Upgrading-Existing-Scripts
# 

# Note: Start Cytoscape
# ------------------------------------------------------------------------
cytoscapePing ()
```

    ## [1] "You are connected to Cytoscape!"

``` r
cytoscapeVersionInfo ()
```

    ##       apiVersion cytoscapeVersion 
    ##             "v1"          "3.7.2"

## Functions for creating and manipulating networks using igraph and RCy3

    ## [1] "RCy3 2.N.N Functions Loaded!"

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
