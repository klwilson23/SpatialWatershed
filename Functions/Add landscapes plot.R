add2plot <- function()
{
  par(xpd=TRUE,usr=c(-1.572561,1.572561,-1.080000,1.080000),mar=c(1.75,1.25,0.75,1.25))
  nodeScalar <- 30
  network <- makeNetworks("linear",Npatches=Npatches,patchDist=patchDist)
  plot(network$landscape,vertex.color="grey",layout=cbind(seq(-1,1,length.out = gorder(network$landscape)),seq(1,-1,length.out = gorder(network$landscape))),vertex.size=network$node.size*nodeScalar,xlim=c(-1,1),ylim=c(-1,1),rescale=FALSE,vertex.label=NA)
  
  nodeScalar <- 30
  network <- makeNetworks("dendritic",Npatches=Npatches,patchDist=patchDist)
  plot(network$landscape,vertex.color="grey",layout=layout_as_tree(network$landscape,root=V(network$landscape)[1]),vertex.size=network$node.size*nodeScalar,vertex.label=NA)
  
  network <- makeNetworks("star",Npatches=Npatches,patchDist=patchDist)
  plot(network$landscape,vertex.color="grey",layout=layout.reingold.tilford(network$landscape,circular=T),vertex.size=network$node.size*nodeScalar,vertex.label=NA)
  
  network <- makeNetworks("complex",Npatches=Npatches,patchDist=patchDist)
  spatialLayout <- matrix(MORELETTERS(1:Npatches),nrow=sqrt(Npatches),ncol=sqrt(Npatches),byrow=TRUE)
  spatialLayout <- t(sapply(1:Npatches,function(x){which(spatialLayout==LETTERS[x],arr.ind=TRUE)}))
  spatialLayout <- spatialLayout[rank(attr(V(network$landscape),"names")),]
  plot(network$landscape,vertex.color="grey",layout=spatialLayout,vertex.size=network$node.size*nodeScalar,vertex.label=NA)
}

add2plotv2 <- function(x1=1,x2=1,x3=1,x4=1)
{
  par(xpd=TRUE,usr=c(-1.572561,1.572561,-1.080000,1.080000),mar=c(x1,x2,x3,x4))
  nodeScalar <- 30
  network <- makeNetworks("linear",Npatches=Npatches,patchDist=patchDist)
  plot(network$landscape,vertex.color="grey",layout=cbind(seq(-1,1,length.out = gorder(network$landscape)),seq(1,-1,length.out = gorder(network$landscape))),vertex.size=network$node.size*nodeScalar,xlim=c(-1,1),ylim=c(-1,1),rescale=FALSE,vertex.label=NA)
  
  nodeScalar <- 30
  network <- makeNetworks("dendritic",Npatches=Npatches,patchDist=patchDist)
  plot(network$landscape,vertex.color="grey",layout=layout_as_tree(network$landscape,root=V(network$landscape)[1]),vertex.size=network$node.size*nodeScalar,vertex.label=NA)
  
  network <- makeNetworks("star",Npatches=Npatches,patchDist=patchDist)
  plot(network$landscape,vertex.color="grey",layout=layout.reingold.tilford(network$landscape,circular=T),vertex.size=network$node.size*nodeScalar,vertex.label=NA)
  
  network <- makeNetworks("complex",Npatches=Npatches,patchDist=patchDist)
  spatialLayout <- matrix(MORELETTERS(1:Npatches),nrow=sqrt(Npatches),ncol=sqrt(Npatches),byrow=TRUE)
  spatialLayout <- t(sapply(1:Npatches,function(x){which(spatialLayout==LETTERS[x],arr.ind=TRUE)}))
  spatialLayout <- spatialLayout[rank(attr(V(network$landscape),"names")),]
  plot(network$landscape,vertex.color="grey",layout=spatialLayout,vertex.size=network$node.size*nodeScalar,vertex.label=NA)
}
