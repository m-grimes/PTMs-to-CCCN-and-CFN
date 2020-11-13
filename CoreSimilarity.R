
# A Degeneracy Framework for Graph Similarity
# by Giannis Nikolentzos, Polykarpos Meladianos, Stratis Limnios and Michalis Vazirgiannis
# https://www.ijcai.org/Proceedings/2018/0360.pdf

# The article lists three kernels, and the best performant is the 
# shortest path kernel (SP) [Borgwardt and Kriegel, 2005]
# which counts pairs of shortest paths in two graphs having the same 
# source and sink labels and identical path lengths on both graphs.
# If a pair of vertices a and b are disconnected on both graphs (i.e., length(a,b)=\infty), 
# we do not consider the pair to have the same length on both graphs.

# input: two igraph objects, and a verbose flag to output info about computations
# output: a similarity value in [0,1]


coreSimilarity<-function(g1, g2,verbose=T){
  require(igraph)
  if(vcount(g1)==0) return(0);
  core1<-graph.coreness(g1);
  if(vcount(g2)==0) return(0);
  core2<-graph.coreness(g2);
  if(verbose){
    message("Graphs have max cores: ",max(core1)," and ",max(core2))
  }
  maxCore = max(core1,core2)
  equiSum <- 0
  pairCount<- 0
  pairEdges <-0
  for( core in 1:maxCore){
    simCore=0
    if(core %in% core1& core %in% core2){
      c1Names = (which(core1==core))
      c2Names = (which(core2==core))
      c1Andc2 = (intersect(names(c1Names),names(c2Names)))
      commonPairs <- length(c1Andc2)
      pairCount = pairCount+commonPairs
      if(commonPairs>1){
        g1Dist = distances(g1, v=c1Andc2,to=c1Andc2,mode="all")
        g2Dist =  distances(g2, v=c1Andc2,to=c1Andc2,mode="all")
        
        for(v1 in 1:(commonPairs-1)){
          for(v2 in (v1+1):commonPairs){
            ver1 = c1Andc2[v1]
            ver2 = c1Andc2[v2]
            dist1 = g1Dist[ver1,ver2]
            dist2 = g2Dist[ver1,ver2]
            if(is.finite(dist1)){
              pairEdges=pairEdges+1
              if(is.finite(dist2)){
                if(dist1==dist2){
                  simCore = simCore+1
                }
              }
            }
          }
        }
        equiSum = equiSum+simCore
        if(verbose){
        message(core," Core: ",commonPairs," common vertices and ",simCore, " equidistant vertex pairs.")
        }
      }
    }
  }
  if(pairEdges==0) resVal=0 
  else resVal=equiSum/(pairEdges)
  if(verbose){
    message("On all cores, ",pairCount," vertices were common in the two graphs.")
    message("A total of ",equiSum," equidistant vertex pairs were found, out of ",pairEdges, " possible.")
    message("Sim = ",equiSum,"/",pairEdges," = ",resVal)
  }
  return(resVal)
}

# example usage
#g1 <- graph_from_literal( Alice-Bob-Cecil-Alice, Daniel-Cecil-Eugene,Cecil-Gordon )
#g2 <- graph_from_literal( Alice-Cuneyt-Cecil-Alice, Daniel-Cecil-Eugene,Cecil-Cemil )
#g3 <- graph_from_literal( Alice-Cuneyt)
#sim<-coreSimilarity(g2,g1,verbose=T)



