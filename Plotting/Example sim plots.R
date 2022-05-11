source("Functions/Linear network.R")
source("Functions/Dispersal function.R")
source("Functions/patch_variance.R")
source("Functions/local disturbance.R")
source("Functions/some functions.R")
source("Functions/Metapop function.R")
source("Functions/popDynFn.R")
source("Functions/Add landscapes plot.R")

matLayout <- matrix(0,ncol=18,nrow=25,byrow=T)
matLayout[1:12,1:9] <- 1
matLayout[1:4,10:12] <- 2
matLayout[1:4,13:15] <- 3
matLayout[5:8,10:12] <- 4
matLayout[5:8,13:15] <- 5
matLayout[9:12,10:12] <- 6
matLayout[9:12,13:15] <- 7
matLayout[7:18,16:17] <- 8

matLayout[14:25,1:9] <- 9
matLayout[14:17,10:12] <- 10
matLayout[14:17,13:15] <- 11
matLayout[18:21,10:12] <- 12
matLayout[18:21,13:15] <- 13
matLayout[22:25,10:12] <- 14
matLayout[22:25,13:15] <- 15
seed_no <- 535
#jpeg("Figures/example landscape results.jpeg",units="in",height=5.5,width=6,res=800)
layout(matLayout)
set.seed(seed_no)
Nburnin<-50
linearSim <- metaPop(Npatches=16,networkType="linear",patchDistance=1,Nburnin=50,NyrsPost=50,omega=0.01,m=100,alpha=2,metaK=1600,cv=1e-1,DistScenario="random_patch",magnitude_of_decline=0.9,lagTime=1,prodType="Beverton-Holt",rho.time=1e-5,rho.dist=1e5,compensationLag=25,dataWeighting=0.1,alphaVariable=FALSE,kVariable=FALSE,spatialPlots=FALSE)
metaK <- 1600
spatialRecoveryPlotv2(textSize=1,linearSim$popDyn,linearSim$MetaPop,linearSim$k_p,Nlevels=10,linearSim$recovery,50,50+50,nodeScalar=35,linearSim$network,networkType="linear",Npatches=Npatches,panel_text=c("(a)","(b)"),year_seq=(Nburnin+seq(0,50,by=10)))
Nlevels <- 10
textSize <- 1
par(mar=c(1,1,1,1))
plot(NA,xlim=c(-1,1),ylim=c(-1,1),xaxt="n",yaxt="n",xlab="",ylab="",frame.plot=FALSE,xpd=NA)
colorlegend(rev(colfunc(Nlevels)),zlevels=Nlevels,zlim=c(0,1),dz=0.1,posx=c(0.35,0.5),posy=c(0.1,0.95),digit=2)
title(main=expression('N'[t]*'/K'),line=0.5,font.main=1,cex.main=1.5*textSize,xpd=NA)
set.seed(seed_no)
dendriticSim <- metaPop(Npatches=16,networkType="dendritic",patchDistance=1,Nburnin=50,NyrsPost=50,omega=0.01,m=100,alpha=2,metaK=1600,cv=1e-1,DistScenario="random_patch",magnitude_of_decline=0.9,lagTime=1,prodType="Beverton-Holt",rho.time=1e-5,rho.dist=1e5,compensationLag=25,dataWeighting=0.1,alphaVariable=FALSE,kVariable=FALSE,spatialPlots=FALSE)
spatialRecoveryPlotv2(textSize=1,dendriticSim$popDyn,dendriticSim$MetaPop,dendriticSim$k_p,Nlevels=10,dendriticSim$recovery,50,50+50,nodeScalar=35,dendriticSim$network,networkType="dendritic",Npatches=Npatches,panel_text=c("(c)","(d)"),year_seq=(Nburnin+seq(0,50,by=10)))
#dev.off()
apply(dendriticSim$popDyn,c(3,2),function(x){(sd(x)/sqrt(1-1e-5^2))/mean(x)})