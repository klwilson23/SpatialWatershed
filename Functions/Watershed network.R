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
NlocalPops <- 4 # how many populations in metapopulation?
Nbranches <- 1 # how many branches does each local pop have?
Nforks <- 4 # how many forks does each branch have?

diameter <- 15
colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))

# declare coordinates of patches for: (1) main branches, (2) secondary branhes, (3) tertiary branches
# add more levels if desired.
edges1 <- sapply(1:NmetaPatch,function(y){sapply(LETTERS[1:NlocalPops],function(x){c(y,x)})})
edges2 <- sapply(LETTERS[1:NlocalPops],function(y){sapply(paste(y,1:Nbranches,sep=""),function(x){c(y,x)})})
edges3 <- sapply(as.vector(unique(edges2[-c(1,3),])),function(y){sapply(paste(y,letters[1:Nforks],sep=""),function(x)c(y,x))})
riverLandscape <- graph(edges=as.character(c(edges1,edges2,edges3)),directed=F)

# set distance between patches for: (1) main branches, (2) secondary branhes, (3) tertiary branches
# add more levels if desired.
dist1 <- rep(diameter,NlocalPops)
dist2 <- rep(diameter*1,(NlocalPops*Nbranches))
dist3 <- rep(diameter*5,(NlocalPops*Nbranches*Nforks))
E(riverLandscape)$weight <- 1/c(dist1,dist2,dist3)

node.size<-setNames(c(2,rep(0.8,NlocalPops),rep(0.8,NlocalPops*Nbranches),rep(0.8,NlocalPops*Nbranches*Nforks)),V(riverLandscape)$names)

V(riverLandscape)$color[c(1,as.vector(t(sapply(LETTERS[1:NlocalPops],function(x){grep(x,V(riverLandscape)$name)}))))] <- c("grey50",as.vector(t(sapply(LETTERS[1:NlocalPops],function(x){rep(colfunc(NlocalPops)[which(x==LETTERS[1:NlocalPops])],length(grep(x,V(riverLandscape)$name)))}))))

#layout(1)
#par(mar=c(5,4,1,1))
#tiff("watershed.tiff",compression="lzw",units="in",height=8,width=8,res=800)
#plot(riverLandscape,col="dodgerblue",layout=layout.auto(riverLandscape),vertex.size=node.size*15)
#dev.off()

distance_table(riverLandscape)
mean_distance(riverLandscape)
distances(riverLandscape,v=V(riverLandscape),to=V(riverLandscape))

