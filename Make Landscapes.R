source("some functions.R")
source("make networks.R")
patchDist <- 1
Npatches <- 16
nodeScalar <- 17
tiff("networks.tiff",compression="lzw",units="in",height=8,width=8,res=800)
layout(matrix(1:4,nrow=2,ncol=2,byrow=T))
par(mar=c(1,1,1,1),oma=c(0,0,0,0))

network <- makeNetworks("linear",Npatches=Npatches,patchDist=patchDist)
plot(network$landscape,col="dodgerblue",layout=cbind(0,seq(1,-1,length.out = gorder(network$landscape))),vertex.size=network$node.size*nodeScalar,xlim=c(-1,1),ylim=c(-1,1),rescale=FALSE)
title(main="Linear",line=0,font.main=1)

network <- makeNetworks("dendritic",Npatches=Npatches,patchDist=patchDist)
plot(network$landscape,col="dodgerblue",layout=layout_as_tree(network$landscape,root=V(network$landscape)[1]),vertex.size=network$node.size*nodeScalar)
title(main="Dendritic",line=0,font.main=1)

network <- makeNetworks("star",Npatches=Npatches,patchDist=patchDist)
plot(network$landscape,col="dodgerblue",layout=layout.reingold.tilford(network$landscape,circular=T),vertex.size=network$node.size*nodeScalar)
title(main="Star",line=0,font.main=1)

network <- makeNetworks("complex",Npatches=Npatches,patchDist=patchDist)
plot(network$landscape,col="dodgerblue",layout=layout_nicely(network$landscape,dim=2),vertex.size=network$node.size*nodeScalar)
title(main="Complex",line=0,font.main=1)
dev.off()

jpeg("networks.jpeg",units="in",height=8,width=8,res=800)
layout(matrix(1:4,nrow=2,ncol=2,byrow=T))
par(mar=c(1,1,1,1),oma=c(0,0,0,0))

network <- makeNetworks("linear",Npatches=Npatches,patchDist=patchDist)
plot(network$landscape,col="dodgerblue",layout=cbind(0,seq(1,-1,length.out = gorder(network$landscape))),vertex.size=network$node.size*nodeScalar,xlim=c(-1,1),ylim=c(-1,1),rescale=FALSE)
title(main="Linear",line=0,font.main=1)

network <- makeNetworks("dendritic",Npatches=Npatches,patchDist=patchDist)
plot(network$landscape,col="dodgerblue",layout=layout_as_tree(network$landscape,root=V(network$landscape)[1]),vertex.size=network$node.size*nodeScalar)
title(main="Dendritic",line=0,font.main=1)

network <- makeNetworks("star",Npatches=Npatches,patchDist=patchDist)
plot(network$landscape,col="dodgerblue",layout=layout.reingold.tilford(network$landscape,circular=T),vertex.size=network$node.size*nodeScalar)
title(main="Star",line=0,font.main=1)

network <- makeNetworks("complex",Npatches=Npatches,patchDist=patchDist)
plot(network$landscape,col="dodgerblue",layout=layout_nicely(network$landscape,dim=2),vertex.size=network$node.size*nodeScalar)
title(main="Complex",line=0,font.main=1)
dev.off()
