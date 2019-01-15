# Author: Kyle L Wilson
# Institution: University of Calgary
# Date: July 26, 2017
# Create simulation social-ecological model based on Hunt et al. 2011
# Date: April 10, 2018
# create 'makeLandscape' function that generates spatial features of the landscape SES model

library(igraph)
library(boot)
library(timeSeries)
library(shape)

source("some functions.R")

NmetaPatch <- 1 # how many metapopulation patches?
NlocalPops <- 24 # how many populations in metapopulation?

MORELETTERS <- extend(LETTERS)

patches <- matrix(c(MORELETTERS(1:(NlocalPops/2)),"1",MORELETTERS((NlocalPops/2+1):NlocalPops)),nrow=5,ncol=5,byrow=T)


diameter <- 15
colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))

# declare coordinates of patches for: level (1) branches, level (2) branhes, level (3) branches, and level (4) branches
# add more levels if desired.

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

# set distance between patches for: (1) main branches, (2) secondary branhes, (3) tertiary branches
# add more levels if desired.
dist2 <- rep(diameter/2,length(edgesTot))
E(complexLandscape)$weight <- 1/dist2

E(complexLandscape)$weight[grep("1",attr(E(complexLandscape),"vnames"))] <- 1/diameter

node.size<-setNames(c(rep(0.8,NlocalPops/2),2,rep(0.8,NlocalPops/2)),V(complexLandscape)$names)

V(complexLandscape)$color[c(which(patches=="1"),as.vector(t(sapply(patches[-which(patches=="1")],function(x){match(x,V(complexLandscape)$name)}))))] <- c("grey50",as.vector(t(sapply(MORELETTERS(1:NlocalPops),function(x){rep(colfunc(NlocalPops)[which(x==MORELETTERS(1:NlocalPops))],length(match(x,V(complexLandscape)$name)))}))))

#layout(1)
#par(mar=c(5,4,1,1))
#tiff("watershed.tiff",compression="lzw",units="in",height=8,width=8,res=800)
#plot(complexLandscape,col="dodgerblue",layout=layout.auto(complexLandscape),vertex.size=node.size*15)
#dev.off()

distance_table(complexLandscape)
mean_distance(complexLandscape)
distances(complexLandscape,v=V(complexLandscape),to=V(complexLandscape))

