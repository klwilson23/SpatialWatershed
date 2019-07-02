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

library(mvtnorm)

source("Functions/Linear network.R")
source("Functions/Dispersal function.R")
source("Functions/patch_variance.R")
source("Functions/local disturbance.R")
source("Functions/some functions.R")
source("Functions/Metapop function.R")
source("Functions/popDynFn.R")
source("Functions/Add landscapes plot.R")

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
dist_colours <- c("blue","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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
dist_colours <- c("blue","dodgerblue","orange")
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