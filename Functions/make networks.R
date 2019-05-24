# Author: Kyle L Wilson
# Institution: Simon Fraser University
# Date: March 29, 2019
# make general function to create landscape and return distance matrix

library(igraph)
library(shape)
source("Functions/some functions.R")

makeNetworks <- function(network,Npatches,patchDist)
{
  MORELETTERS <- extend(LETTERS)
  colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))
  diameter <- patchDist
  
  if(network=="linear")
  {
    edges1 <- as.vector(t(cbind(MORELETTERS(1:(Npatches-1)),MORELETTERS(1:Npatches)[-1])))
    linearLandscape <- graph(edges=as.character(edges1),directed=F)
    # set distance between patches
    # add more levels if desired.
    dist1 <- rep(diameter,length(as.character(edges1))/2)
    E(linearLandscape)$weight <- 1/dist1
    node.size<-setNames(rep(0.8,Npatches),V(linearLandscape)$names)
    V(linearLandscape)$color[as.vector(t(sapply(MORELETTERS(1:Npatches),function(x){match(x,V(linearLandscape)$name)})))] <- as.vector(t(sapply(MORELETTERS(1:Npatches),function(x){rep(colfunc(Npatches)[which(x==MORELETTERS(1:Npatches))],length(match(x,V(linearLandscape)$name)))})))
    dist_matrix <- distances(linearLandscape,v=V(linearLandscape),to=V(linearLandscape))
    landscape <- linearLandscape
  }
  
  if(network=="star")
  {
    NcentralPatch <- 1 # how many central patches?
    NlocalPops <- Npatches-NcentralPatch # how many populations in metapopulation? should match up in total patches as the Watershed Network
    # declare arrangement of patches around metapopulation
    edges1 <- sapply(MORELETTERS(1:NcentralPatch),function(y){sapply(MORELETTERS((NcentralPatch+1):Npatches),function(x){c(y,x)})})
    starLandscape <- graph(edges=as.character(edges1),directed=F)
    
    # set distance between patches around metapopulation
    # add more levels if desired.
    dist1 <- rep(diameter,NlocalPops)
    E(starLandscape)$weight <- 1/dist1
    node.size<-setNames(rep(0.8,Npatches),V(starLandscape)$names)
    V(starLandscape)$color[as.vector(t(sapply(MORELETTERS(1:Npatches),function(x){match(x,V(starLandscape)$name)})))] <- as.vector(t(sapply(MORELETTERS(1:Npatches),function(x){rep(colfunc(Npatches)[which(x==MORELETTERS(1:Npatches))],length(match(x,V(starLandscape)$name)))})))
    dist_matrix <- distances(starLandscape,v=V(starLandscape),to=V(starLandscape))
    landscape <- starLandscape
  }
  
  if(network=="dendritic")
  {
    NmetaPatch <- 1 # how many metapopulation patches?
    NlocalPops <- 1 # how many populations in metapopulation?
    Nbranches <- 2 # how many branches does each local pop have?
    Nforks <- 2 # how many forks does each branch have?
    Ntails <- 2 # how many tails does each fork have?
    
    # declare coordinates of patches for: (1) main branches, (2) secondary branhes, (3) tertiary branches
    # add more levels if desired.
    edges1 <- sapply(MORELETTERS(1:NmetaPatch),function(y){sapply(MORELETTERS((NmetaPatch+1):(NmetaPatch+NlocalPops)),function(x){c(y,x)})})
    edges2 <- sapply(MORELETTERS((NmetaPatch+1):(NmetaPatch+NlocalPops)),function(y){sapply(paste(y,1:Nbranches,sep=""),function(x){c(y,x)})})
    edges3 <- sapply(as.vector(unique(edges2[(1:nrow(edges2)%%2==0),])),function(y){sapply(paste(y,letters[1:Nforks],sep=""),function(x)c(y,x))})
    edges4 <- sapply(as.vector(unique(edges3[(1:nrow(edges2)%%2==0),])),function(y){sapply(paste(y,1:Ntails,sep=""),function(x)c(y,x))})
    riverLandscape <- graph(edges=as.character(c(edges1,edges2,edges3,edges4)),directed=FALSE)
    
    # set distance between patches for: (1) main branches, (2) secondary branhes, (3) tertiary branches
    # add more levels if desired.
    dist1 <- rep(diameter,NmetaPatch*NlocalPops)
    dist2 <- rep(diameter,(NlocalPops*Nbranches))
    dist3 <- rep(diameter,(NlocalPops*Nbranches*Nforks))
    dist4 <- rep(diameter,(NlocalPops*Nbranches*Nforks*Ntails))
    E(riverLandscape)$weight <- 1/c(dist1,dist2,dist3,dist4)
    
    node.size<-setNames(rep(0.8,Npatches),V(riverLandscape)$names)
    
    dist_matrix <- distances(riverLandscape,v=V(riverLandscape),to=V(riverLandscape))
    patchRank <- as.numeric(as.factor(rank(nchar(colnames(dist_matrix)))))
    patchRank[2:length(patchRank)] <- 1+patchRank[2:length(patchRank)]
    patchColor <- rank(colnames(dist_matrix))
    #V(riverLandscape)$color <- colfunc(length(unique(patchRank)))[patchRank]
    V(riverLandscape)$color <- colfunc(length(unique(patchColor)))[patchColor]
    V(riverLandscape)$label <- MORELETTERS(1:Npatches)
    landscape <- riverLandscape
  }
  
  if(network=="complex")
  {
    patches <- matrix(MORELETTERS(1:Npatches),nrow=sqrt(Npatches),ncol=sqrt(Npatches),byrow=T)
    
    edge1 <- c()
    for(i in 1:nrow(patches))
    {
      for(j in 1:ncol(patches))
      {
        edgeList <- as.vector(t(patches[max(1,(i-1)):min(nrow(patches),(i+1)),max(1,(j-1)):min(ncol(patches),(j+1))]))[-which(as.vector(t(patches[max(1,(i-1)):min(nrow(patches),(i+1)),max(1,(j-1)):min(ncol(patches),(j+1))]))==patches[i,j])]
        edge1 <- c(edge1,as.vector(sapply(edgeList,function(x){c(patches[i,j],x)})))
      }
    }
    
    edgesTot <- edge1
    complexLandscape <- graph(edges=as.character(edgesTot),directed=F)
    dist1 <- rep(diameter,length(edgesTot)/2)
    E(complexLandscape)$weight <- 1/dist1
    node.size<-setNames(rep(0.8,Npatches),V(complexLandscape)$names)
    V(complexLandscape)$color <- colfunc(Npatches)
    dist_matrix <- distances(complexLandscape,v=V(complexLandscape),to=V(complexLandscape))
    landscape <- complexLandscape
  }
  
  return(list("distanceMatrix"=dist_matrix,"landscape"=landscape,"node.size"=node.size))
}
