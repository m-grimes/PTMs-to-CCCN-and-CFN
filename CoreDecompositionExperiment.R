library(igraph)
library(ggplot2)
library(reshape)
library(RCy3)
library(pheatmap)
library(tidyr)
library(dplyr)
rm(list = ls())
# Learn where data files are stored.
source("code/dataConfig.R")
#load utility functions
source("code/utilityFunctions.R")
# We compute similarity of two pathways by using different gene similarity measures

load(bioPlanetsRdata)

picDir="figures/"
gzall<-simplify(gzalltgene.physical.cfn.merged.g) 
gldgene<-simplify(ldgene.physical.cfn.merged.g)

commonV<-intersect(names(V(gldgene)),names(V(gzall)))
length(commonV)

egzall<-as.data.frame(as_edgelist(gzall, names = TRUE))
egldgene<-as.data.frame(as_edgelist(gldgene, names = TRUE))
commonE<-inner_join(egzall, egldgene)


egzall2<-egzall[egzall$V1%in%commonV&egzall$V2%in%commonV,]
egldgene2<-egldgene[egldgene$V1%in%commonV&egldgene$V2%in%commonV,]
allPossibleEdges<-unique(bind_rows(egzall2,egldgene2))
nrow(egzall2)
nrow(egldgene2)
nrow(allPossibleEdges)
nrow(commonE)

edge_density(gzall, loops=F)
edge_density(gldgene, loops=F)
m1<-motifs(gzall, size = 4)
m2<-motifs(gldgene, size = 4)
m1<-m1/sum(m1,na.rm=TRUE)
m2<-m2/sum(m2,na.rm=TRUE)
cosineSim(m1,m2)
m2
# compare high cores
gzallCoreness <- graph.coreness(gzall);
p1<-ggplot() + theme_bw()+
  geom_histogram(aes(x=gzallCoreness),binwidth = 1,color="white", fill="steelblue2")+
scale_x_continuous(name="core",breaks = c(1,4,6,8,10,12))+
  theme_minimal()+scale_y_continuous(name=("#PTMs"))+
  theme(text = element_text(size=22));p1
ggsave(plot=p1,file="gzallCoreHistogram.png",device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)

threshold = 10
gzall_2<-delete_vertices(gzall, names(gzallCoreness[gzallCoreness<threshold]))
ecount(gzall_2)
vcount(gzall_2)

gldcoreness <- graph.coreness(gldgene);
p2<-ggplot() + theme_bw()+
  geom_histogram(aes(x=gldcoreness),binwidth = 1,color="white", fill="steelblue2")+
  scale_x_continuous(name="core",breaks = c(1,4,6,8,9,10))+
  theme_minimal()+scale_y_continuous(name=("#PTMs"))+
  theme(text = element_text(size=22));p2
ggsave(plot=p2,file="gldgeneCoreHistogram.png",device="png",width=6,height=4,units=c("in"),dpi=1200,path=picDir)
threshold=8
gld_2<-delete_vertices(gldgene, names(gldcoreness[gldcoreness<threshold]))

ecount(gld_2)
vcount(gld_2)
##
# common vertices that appear in high cores
intersect(commonV,names(V(gld_2)))
intersect(commonV,names(V(gzall_2)))

# Check edge weights
# directed edges

deg <- degree(gzall, mode="all")
deg.dist <- degree_distribution(gzall, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, 
      col="orange",main="gzall",xlab="Degree", ylab="Cumulative Frequency")

deg2 <- degree(gldgene, mode="all")
deg2.dist <- degree_distribution(gldgene, cumulative=T, mode="all")
plot( x=0:max(deg2), y=1-deg2.dist, pch=19, cex=1.2, 
      col="orange",main="ldgene",xlab="Degree", ylab="Cumulative Frequency")
message("V count: ",length(V(gzall_2))," ",length(V(gld_2)))
#plot(g1_2, layout=CorenessLayout(g1_2),edge.width=1,edge.arrow.size=0.5, arrow.mode=0,vertex.size=1, vertex.color=colbar[coreness], vertex.frame.color=colbar[coreness], main='Coreness')

## hubs and authorities

hs <- hub_score(gzall)$vector
hs2 <- hub_score(gldgene)$vector

mean_distance(gzall, directed=T)
mean_distance(gldgene, directed=T)

hist(distances(gzall))
mean_distance((gzall))
hist(distances(gldgene))
mean_distance((gldgene))

cfg <- cluster_fast_greedy(simplify(as.undirected(gzall_2)))
plot(cfg, simplify(as.undirected(gzall_2)))

cfg <- cluster_fast_greedy(simplify(as.undirected(gld_2)))
plot(cfg, simplify(as.undirected(gld_2)))

# we manually set the highest cores 
gzallHighCoreNodes<-names(gzallCoreness[gzallCoreness==max(gzallCoreness)])
gldhighCoreNodes<-names(gldcoreness[gldcoreness==max(gldcoreness)])

gzallHighCoreNetwork<-simplify(induced_subgraph(gzall,gzallHighCoreNodes))
cfg <- cluster_fast_greedy(gzallHighCoreNetwork)
plot(cfg, gzallHighCoreNetwork)

gldHighCorenetwork<-simplify(induced_subgraph(gldgene,gldhighCoreNodes))
cfg <- cluster_fast_greedy(gldHighCorenetwork)
plot(cfg, gldHighCorenetwork)

intersect(gzallHighCoreNodes,gldhighCoreNodes)

corematrix <- array(c(0, 0), dim = c(max(gzallCoreness),max(gldcoreness)))
for(i in seq(1,max(gzallCoreness),1)){
  gzallHighCoreNodes<-names(gzallCoreness[gzallCoreness==i])
  for(j in seq(1,max(gldcoreness),1)){
    gldhighCoreNodes<-names(gldcoreness[gldcoreness==j])
    lend=length(intersect(gzallHighCoreNodes,gldhighCoreNodes))
    if(lend==0)lend=NA
    corematrix[i,j]=lend
  }
}
corematrix2<-corematrix
corematrix2[is.na(corematrix2)] <- ""
colnames(corematrix2) <- seq(1,ncol(corematrix),1)
rownames(corematrix2) <- seq(1,nrow(corematrix),1)
png(file="figures/corematrix.png", width = 4, height = 4, units = 'in', res = 2000)
pheatmap(corematrix,border_color = NA,number_format="%.0f",xTitle = "xTest",
         yTitle = "yTest",cluster_rows = F,cluster_cols =F,
         fontsize = 22,legend=F,na_col = "gray",display_numbers = corematrix2)
dev.off()

 

# article pictures from August 2021
gzallHighCoreNodes<-names(gzallCoreness[gzallCoreness>10])

gzallHighCoreNetwork<-simplify(induced_subgraph(gzall,gzallHighCoreNodes))

gldhighCoreNodes<-names(gldcoreness[gldcoreness>6])
gldHighCorenetwork<-simplify(induced_subgraph(gldgene,gldhighCoreNodes))

### common vertices
common<-intersect(gzallHighCoreNodes,gldhighCoreNodes)

gzallHighCoreNetwork<-simplify(induced_subgraph(gzall,common))

V(gzallHighCoreNetwork)$size <-  1
V(gzallHighCoreNetwork)$shape <-"none"  
V(gzallHighCoreNetwork)$label.cex =1
V(gzallHighCoreNetwork)$label.distance =100

cfg <- cluster_fast_greedy(gzallHighCoreNetwork)
cfg <- cluster_louvain(gzallHighCoreNetwork)
cfg <- cluster_edge_betweenness(gzallHighCoreNetwork)
plot(cfg,gzallHighCoreNetwork, margin=c(0,0,0,0))
#png( "gzallHighCoreNetwork.png", 1000, 1000)
plot(cfg,gzallHighCoreNetwork, vertex.size=8, vertex.color = rainbow(10, .8, .8, alpha= .8),
     vertex.label.color = "black", vertex.label.cex = 1, vertex.label.degree = 1,
     edge.arrow.size = 0, edge.arrow.width = 0, edge.color = "gray55")
#dev.off()

gldHighCorenetwork<-simplify(induced_subgraph(gldgene,common))

V(gldHighCorenetwork)$size <-  1
V(gldHighCorenetwork)$shape <-"none"  
V(gldHighCorenetwork)$label.cex =1
V(gldHighCorenetwork)$label.distance =100

cfg <- cluster_fast_greedy(gldHighCorenetwork)
cfg <- cluster_louvain(gldHighCorenetwork)
cfginsta <- cluster_edge_betweenness(gldHighCorenetwork)
plot(cfg,gldHighCorenetwork, margin=c(0,0,0,0))
#png( "gldHighCorenetwork.png", 1000, 1000)
plot(cfg,gldHighCorenetwork, vertex.size=8, vertex.color = rainbow(10, .8, .8, alpha= .8),
     vertex.label.color = "black", vertex.label.cex = 1, vertex.label.degree = 1,
     edge.arrow.size = 0, edge.arrow.width = 0, edge.color = "gray55")
#dev.off()

write_graph(gldHighCorenetwork, "gld.edges", format = "ncol")
write_graph(gzallHighCoreNetwork, "gzall.edges", format = "ncol")


gz.edges <- data.frame(as_edgelist(gzallHighCoreNetwork))
ld.edges <- data.frame(as_edgelist(gldHighCorenetwork))
names(gz.edges) <- c("source", "target")
names(ld.edges) <- c("source", "target")
gz.edges$interaction <- "igraph"
gz.edges$Weight <- 0.2
ld.edges$interaction <- "igraph"
ld.edges$Weight <- 0.2

gz.cf <- data.frame(id=names(V(gzallHighCoreNetwork)))
ld.cf <- data.frame(id=names(V(gldHighCorenetwork)))

tryit.graph <- createNetworkFromDataFrames(gz.cf, gz.edges, "Test Cy from igraph")
tryit.graph <- createNetworkFromDataFrames(ld.cf, ld.edges, "gldgene")
# common edges in the two high core networks 
edges<-intersection(E(gzallHighCoreNetwork),E(gldHighCorenetwork))
nodes<-intersection(V(gzallHighCoreNetwork),V(gldHighCorenetwork))
ed<-data.frame(id=as_ids(edges))

ed=separate(ed, col = id, into = c("source","target"), sep = "\\|")
vd<-data.frame(id=as_ids(nodes))
commongraph <- createNetworkFromDataFrames(vd, ed, "")
# differing edges in the two high core networks

nodes<-intersection(V(gzallHighCoreNetwork),V(gldHighCorenetwork))
gz2<-induced_subgraph(gzallHighCoreNetwork,vids=names(nodes))
edges<-difference(E(gz2),E(gldHighCorenetwork))

ed<-data.frame(id=as_ids(edges))
ed2=separate(ed, col = id, into = c("source","target"), sep = "\\|")
vd<-data.frame(id=as_ids(nodes))
differinggraph <- createNetworkFromDataFrames(vd, ed2, "")

