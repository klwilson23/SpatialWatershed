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

source("Figures/some functions.R")

NmetaPatch <- 1 # how many metapopulation patches?
NlocalPops <- 24 # how many populations in metapopulation? should match up in total patches as the Watershed Network

MORELETTERS <- extend(LETTERS)

diameter <- 15 # distance between patches
colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))

# declare arrangement of patches around metapopulation
# add more levels if desired.
edges1 <- sapply(1:NmetaPatch,function(y){sapply(MORELETTERS(1:NlocalPops),function(x){c(y,x)})})
starLandscape <- graph(edges=as.character(edges1),directed=F)

# set distance between patches around metapopulation
# add more levels if desired.
dist1 <- rep(diameter,NlocalPops)
E(starLandscape)$weight <- 1/dist1

node.size<-setNames(c(2,rep(0.8,NlocalPops)),V(starLandscape)$names)

V(starLandscape)$color[c(1,as.vector(t(sapply(MORELETTERS(1:NlocalPops),function(x){match(x,V(starLandscape)$name)}))))] <- c("grey50",as.vector(t(sapply(MORELETTERS(1:NlocalPops),function(x){rep(colfunc(NlocalPops)[which(x==MORELETTERS(1:NlocalPops))],length(match(x,V(starLandscape)$name)))}))))

#layout(1)
#par(mar=c(5,4,1,1))
#tiff("linear.tiff",compression="lzw",units="in",height=8,width=8,res=800)
#plot(starLandscape,col="dodgerblue",layout=layout.auto(linearLandscape),vertex.size=node.size*15)
#dev.off()

distance_table(starLandscape)
mean_distance(starLandscape)
distances(starLandscape,v=V(starLandscape),to=V(starLandscape))
