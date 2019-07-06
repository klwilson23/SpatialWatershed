results <- readRDS("Simulations/results2019-06-30.rds")
scenarios <- readRDS("Simulations/scenarios2019-06-30.rds")

NaNsims <- which(apply(apply(results[,-(1:8)],2,is.nan),1,any))

results[NaNsims,]

Nboots <- 100
Nlevels <- 11
model <- "Beverton-Holt"
# simulation parameters
Npatches <- 16
patchDist <- 1
Nburnin <- 50
NyrsPost <- 100
Nyears <- Nburnin+NyrsPost
compLag <- 25 # how many years to lag estimates of compensation
dataWeighting <- 0.1 # penalty to past years for data weighting
m <- 100 # distance decay function: penalty of 1 says about 60% of dispersing recruits move to neighbor patches. penalty of 2 says about 90% of dispersing recruits move to neighbor patches
# adult stock-juvenile recruitment traits
alpha <- 2
metaK <- 1600
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.9
# what is the lag time between recruits and spawners
lagTime <- 1
results$StateShift <- ((1-results$recovered)+(1-results$longOcc))/2
results$RecoveryRate <- 1-results$recovery/NyrsPost

library(mvtnorm)
library(vioplot)
source("Functions/Linear network.R")
source("Functions/Dispersal function.R")
source("Functions/patch_variance.R")
source("Functions/local disturbance.R")
source("Functions/some functions.R")
source("Functions/Metapop function.R")
source("Functions/popDynFn.R")
source("Functions/Add landscapes plot.R")

findHull <- function(df,x,y)
{
  df[chull(df[,x],df[,y]),c(x,y)]
}

dist_colours <- c("tomato","dodgerblue","orange")
transparency <- 0.6
plotColours <- ifelse(results$dispersal==0,
                           adjustcolor(dist_colours[results$disturbance],1,offset=c(-0.35,-0.35,-0.35,0)),
                           adjustcolor(dist_colours[results$disturbance],1))
tiff("Figures/Disturbance impacts on recovery regime polygon.tiff",compression="lzw",units="in",height=8,width=8,res=800)
layout(matrix(1:4,nrow=2,ncol=2,byrow=TRUE))
par(mar=c(4,4,0.5,0.5))
plot(results$RecoveryRate,results$longOcc,pch=NA,bg=plotColours,xlab="Recovery rate (yr-1)",ylab="Patch occupancy")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="RecoveryRate",y="longOcc"),col=adjustcolor(dist_colours[z],transparency))}))

plot(results$RecoveryRate,results$longMSY,pch=NA,bg=plotColours,xlab="Recovery rate (yr-1)",ylab="Surplus production")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="RecoveryRate",y="longMSY"),col=adjustcolor(dist_colours[z],transparency))}))

plot(results$longMSY,results$longOcc,pch=NA,bg=plotColours,xlab="Surplus production",ylab="Patch occupancy")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="longMSY",y="longOcc"),col=adjustcolor(dist_colours[z],transparency))}))

plot(results$RecoveryRate,results$recovered,pch=NA,bg=plotColours,xlab="Recovery rate (yr-1)",ylab="Recovered")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="RecoveryRate",y="recovered"),col=adjustcolor(dist_colours[z],transparency))}))

legend("bottomright",c("Uniform","Local, random","Local, extirpation"),pch=22,pt.bg=dist_colours,bty="n",title="Disturbance regime")
dev.off()

boxplot(longCV~stochastic+disturbance,data=results[results$spatial==scenarios$spatial[1] & results$temporal==scenarios$temporal[1],],col=dist_colours[results$disturbance],xlab="Stochasticity",ylab="Recovery rate (yr-1)")

tiff("Figures/Recovery along dispersal network and disturbance clines.tiff",compression="lzw",units="in",height=6,width=6,res=800)

matLayout <- matrix(0,nrow=16,ncol=16)
matLayout[1:8,1:8] <- 1
matLayout[1:8,9:16] <- 2
matLayout[9:16,1:8] <- 3
matLayout[9:16,9:16] <- 4
matLayout[1:3,6:8] <- 5
matLayout[1:3,14:16] <- 6
matLayout[9:11,6:8] <- 7
matLayout[9:11,14:16] <- 8

layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery (generations)",xlab="Dispersal rate",ylim=1.5*range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    #lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=subResults[subResults$disturbance==scenarios$disturbance[z],],x="dispersal",y="recovery"),col=adjustcolor(dist_colours[z],transparency))}))
  abline(h=25,lwd=2,col="black",xpd=FALSE)
}
add2plot()
dev.off()

layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i],]
  vioplot(recovery~alpha,data=subResults,ylim=1.7*range(c(0,results$recovery)),drawRect=FALSE,col=c("orange","dodgerblue"))
  boxplot(recovery~alpha,data=subResults,ylim=1.7*range(c(0,results$recovery)),add=TRUE,col="grey30",boxwex=0.25,xaxt="n",yaxt="n",border="white",outpch=21,outcol="white",whisklty="solid",staplelwd=NA,outcex=1.5,outbg="grey30")
  abline(h=25,lwd=2,col="black",xpd=FALSE)
  mtext("Local productivity",1,line=2.5,cex=0.7)
  mtext("Recovery (generations)",2,line=2.5,cex=0.7)
}
add2plot()

tiff("Figures/Recovery along variable local productivities.tiff",compression="lzw",units="in",height=6,width=6,res=800)
layout(1)
par(mar=c(5,4,1,1))
vioplot(recovery~alpha,data=results,ylim=range(c(0,results$recovery)),drawRect=FALSE,col=c("orange","dodgerblue"))
boxplot(recovery~alpha,data=results,ylim=range(c(0,results$recovery)),add=TRUE,col="grey30",boxwex=0.25,xaxt="n",yaxt="n",border="white",outpch=21,outcol="white",whisklty="solid",staplelwd=NA,outcex=1.5,outbg="grey30")
abline(h=25,lwd=2,col="black",xpd=FALSE)
mtext("Local productivity",1,line=2.5,cex=0.9)
mtext("Recovery (generations)",2,line=2.5,cex=0.9)
dev.off()

layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery (generations)",xlab="Dispersal rate",ylim=1.5*range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery (generations)",xlab="Dispersal rate",ylim=1.5*range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[3] & 
                          results$spatial==scenarios$spatial[3] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery (generations)",xlab="Dispersal rate",ylim=1.5*range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY patterns
# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[2] & 
                          results$spatial==scenarios$spatial[2] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# Capacity
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$med_capacity,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Carrying capacity",xlab="Dispersal rate",ylim=1.5*range(c(0,results$med_capacity)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$med_capacity,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# Compensation
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$med_compensation,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Compensation",xlab="Dispersal rate",ylim=1.5*range(c(0,results$med_compensation)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$med_compensation,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# Occupancy
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[2] & 
                          results$spatial==scenarios$spatial[2] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medOcc,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Patch occupancy",xlab="Dispersal rate",ylim=1.70*range(c(0,results$medOcc)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medOcc,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()