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

rm(list=ls(all=T))

Corner_text <- function(text, location="topright",...) #function to write text to the corner of plots
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}

get_beta <- function(mean,cv) #function that returns the alpha and beta shape parameters of a beta distribution, based on the mean and variation of a given beta distribution
{
  sd <- mean*cv
  alpha <- -((mean*(mean^2+sd^2-mean))/sd^2)
  beta <- alpha/mean-alpha
  return(list(alpha=alpha,beta=beta))
}

layout_in_circles <- function(g, group=1) {
  layout <- lapply(split(V(g), group), function(x) {
    layout_in_circle(induced_subgraph(g,x))
  })
  layout <- Map(`*`, layout, seq_along(layout))
  x <- matrix(0, nrow=vcount(g), ncol=2)
  split(x, group) <- layout
  x
}

NmetaPatch <- 1 # how many towns: 1, 2, 3, or 4
NlocalPops <- 4 # how many lakes in the population
Nbranches <- 2 # how many branches does each local pop have?
Nforks <- 2 # how many forks does each branch have?

migrating <- F # do fish populations migrate between populations or not?
diameter <- 15
maxDist <- 10 # maximum distance of the landscape
colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))
Nyears <- 50 # how many years of the simulation
Nboot <- 10 # how many iterations to run

edges1 <- sapply(1:NmetaPatch,function(y){sapply(LETTERS[1:NlocalPops],function(x){c(y,x)})})
edges2 <- sapply(LETTERS[1:NlocalPops],function(y){sapply(paste(y,1:Nbranches,sep=""),function(x){c(y,x)})})
edges3 <- sapply(as.vector(unique(edges2[-c(1,3),])),function(y){sapply(paste(y,letters[1:Nforks],sep=""),function(x)c(y,x))})
g1 <- graph(edges=as.character(c(edges1,edges2,edges3)),directed=F)

dist1 <- rep(diameter,NlocalPops)
dist2 <- rep(diameter/2,(NlocalPops*Nbranches))
dist3 <- rep(diameter/3,(NlocalPops*Nbranches*Nforks))
E(g1)$weight <- c(dist1,dist2,dist3)

node.size<-setNames(c(2,rep(0.8,NlocalPops),rep(0.8,NlocalPops*Nbranches),rep(0.8,NlocalPops*Nbranches*Nforks)),V(g1)$names)

colfunc(1)

V(g1)$color[c(1,as.vector(t(sapply(LETTERS[1:NlocalPops],function(x){grep(x,V(g1)$name)}))))] <- c("grey50",as.vector(t(sapply(LETTERS[1:NlocalPops],function(x){rep(colfunc(NlocalPops)[which(x==LETTERS[1:NlocalPops])],length(grep(x,V(g1)$name)))}))))

layout(1)
par(mar=c(5,4,1,1))
tiff("watershed.tiff",compression="lzw",units="in",height=8,width=8,res=800)
plot(g1,col="dodgerblue",layout=layout.reingold.tilford(g1,circular=T),vertex.size=node.size*15)
dev.off()

distance_table(g1)

mean_distance(g1)

distances(g1,v=V(g1),to=V(g1))

