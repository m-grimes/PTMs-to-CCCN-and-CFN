# Workaround function for RCy3 glitches
delete.bad.networks <- function(){
  pingtest=cytoscapePing()
  if (!grepl("connected", pingtest)) {
    print("You are NOT conncected to Cytoscape!")
    break } else {
      networks <- getNetworkList()
      for (i in 1:length(networks)) {
        edgeTable <- getTableColumns("edge", c("name", "shared name"), network=networks[i])
        if(!identical(edgeTable[,1], edgeTable[,2])) {
          print(paste('Network', networks[i], "is bad."))
          deleteNetwork(networks[i])} else {
            print (paste("Network", networks[i], "passes edge test."))
          } 
      }
    }
}

# Workaround for edge mapping problem in Cytoscape
FixEdgeDprops.RCy32 <- function() {
  setEdgeLineWidthDefault (3)
  setEdgeColorDefault ( "#FFFFFF")  # white
  setEdgeSelectionColorDefault (col2hex("chartreuse"))
  edgecolors <- col2hex(c("red", "red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2", "darkorange1"))
  edgecolorsplus <- col2hex(c("deeppink", "red", "red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2", "orangered2", "darkorange1"))
  #  red; turquois; green; magenta; blue; violet; green;  bluegreen; black; gray; turquoiseblue; orange 
  edgeTypes <- c("pp", "PHOSPHORYLATION", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of",  "controls-state-change-of", "Physical interactions", "BioPlex", "in-complex-with",  'experiments',  'database',   "Pathway", "Predicted", "Genetic interactions", "correlation", "negative correlation", "positive correlation",  'combined_score', "merged" , "intersect", "peptide", 'homology', "Shared protein domains", "ACETYLATION")
  # 22 edgeTypes            
  myarrows <- c ('Arrow', 'Arrow', 'Arrow','Arrow', 'Arrow', "Arrow", 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', "Arrow")
  setEdgeTargetArrowMapping( 'shared interaction', edgeTypes, myarrows, default.shape='None')  
  matchArrowColorToEdge('TRUE')
  setEdgeColorMapping( 'shared interaction', edgeTypes, edgecolors, 'd', default.color="#FFFFFF")
  # A grey background helps highlight some of the edges
  setBackgroundColorDefault("#949494") # grey 58
  edgevalues <- getTableColumns('edge')
  if (length(edgevalues[grep("pp", edgevalues$'shared interaction'), 1])>0) {
    setEdgeColorBypass(edgevalues[grep("pp", edgevalues$'shared interaction'), 1], col2hex("red"))}
  if (length(edgevalues[grep("PHOSPHORYLATION", edgevalues$'shared interaction'), 1])>0) {
    setEdgeColorBypass(edgevalues[grep("PHOSPHORYLATION", edgevalues$'shared interaction'), 1], col2hex("red"))}
  if (length(edgevalues[grep("phosphorylation", edgevalues$'shared interaction'), 1])>0) {
    setEdgeColorBypass(edgevalues[grep("phosphorylation", edgevalues$'shared interaction'), 1], col2hex("red"))}
  if (length(edgevalues[grep("ACETYLATION", edgevalues$'shared interaction'), 1])>0) {
    setEdgeColorBypass(edgevalues[grep("ACETYLATION", edgevalues$'shared interaction'), 1], col2hex("darkorange1"))}
}

# Function to set edge appearance - revised to test for mapping problem
setCorrEdgeAppearance <- function(edgefile) {
  setEdgeLineWidthDefault (3)
  setEdgeColorDefault ( "#FFFFFF")  # white
  edgevalues <- getTableColumns('edge',c('Weight'))
  edgevalues['Weight']<-abs(edgevalues['Weight'])
  edgevalues['Weight']<-lapply(edgevalues['Weight'], function(x) log2(x * 10) + 2)
  #setEdgeLineWidthBypass(edgevalues[['name']], edgevalues[['Weight']])
  names(edgevalues)<-c('Width')
  loadTableData(edgevalues, table = 'edge', table.key.column = 'SUID')
  setEdgeLineWidthMapping('Width', mapping.type = 'passthrough', style.name = 'default')
  setEdgeSelectionColorDefault ( "#FF69B4")  # hotpink
  edgecolors <- col2hex(c("red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2"))
  edgecolorsplus <- col2hex(c("deeppink", "red", "red", "magenta", "violet", "purple",  "green", "green2", "green3",  "aquamarine2", "cyan", "turquoise2", "cyan2", "lightseagreen", "gold",  "blue", "yellow", "slategrey", "darkslategrey", "grey", "black", "orange", "orange2", "orangered2"))
  #  red; turquois; green; magenta; blue; violet; green;  bluegreen; black; gray; turquoiseblue; orange 
  edgeTypes <- c("pp", "controls-phosphorylation-of", "controls-expression-of", "controls-transport-of",  "controls-state-change-of", "Physical interactions", "BioPlex", "in-complex-with",  'experiments',  'database',   "Pathway", "Predicted", "Genetic interactions", "correlation", "negative correlation", "positive correlation",  'combined_score', "merged" , "intersect", "peptide", 'homology', "Shared protein domains") 
  # 22 edgeTypes            
  myarrows <- c ('Arrow', 'Arrow', 'Arrow', 'Arrow', "Arrow", 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None')
  setEdgeTargetArrowMapping( 'interaction', edgeTypes, myarrows, default.shape='None')  
  matchArrowColorToEdge('TRUE')
  # work around  misaligned mapping
  test <- getTableColumns('edge',c('interaction', "shared interaction", "shared name"))
  if (identical(test$interaction, test$`shared interaction`)) {
    setEdgeColorMapping( 'interaction', edgeTypes, edgecolors, 'd', default.color="#FFFFFF") 
    edgeDprops.RCy32()
  } else {
    setEdgeColorMapping( 'shared interaction', edgeTypes, edgecolors, 'd', default.color="#FFFFFF") 
    FixEdgeDprops.RCy32() }
} 


# GitHup token update: see: https://stackoverflow.com/questions/66065099/how-to-update-github-authentification-token-on-rstudio-to-match-the-new-policy
install.packages("gitcreds")
library(gitcreds)
gitcreds_set()
