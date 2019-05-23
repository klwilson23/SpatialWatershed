results <- readRDS("Simulations/results2019-05-19_fixed.rds")
scenarios <- readRDS("Simulations/scenarios2019-05-19.rds")

NaNsims <- which(apply(apply(results[,-(1:8)],2,is.nan),1,any))

results[NaNsims,]

Nboots <- 50
Nlevels <- 11
model <- "Beverton-Holt"
# simulation parameters
Npatches <- 16
patchDist <- 1
Nburnin <- 50
NyrsPost <- 50
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

source("Linear network.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("Metapop function.R")
source("popDynFn.R")

matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery (generations)",xlab="Dispersal rate",ylim=range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")

matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery (generations)",xlab="Dispersal rate",ylim=range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")

matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery (generations)",xlab="Dispersal rate",ylim=range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")

matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[3] & 
                          results$spatial==scenarios$spatial[3] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery (generations)",xlab="Dispersal rate",ylim=range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")

# MSY patterns

# MSY
matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Maximum sustainable loss",xlab="Dispersal rate",ylim=range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")


# MSY
matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Maximum sustainable loss",xlab="Dispersal rate",ylim=range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")

# MSY
matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[3] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Maximum sustainable loss",xlab="Dispersal rate",ylim=range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")

# MSY
matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$longMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Maximum sustainable loss",xlab="Dispersal rate",ylim=range(c(0,results$longMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$longMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")

# MSY
matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[2] & 
                          results$spatial==scenarios$spatial[2] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medCV,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Medium term Variance",xlab="Dispersal rate",ylim=range(c(0,pmax(0,results$medCV,na.rm=TRUE))))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medCV,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")


# Capacity
matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$med_capacity,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Carrying capacity",xlab="Dispersal rate",ylim=range(c(0,pmax(0,results$med_capacity,na.rm=TRUE))))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$med_capacity,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")

#
# MSY
matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange")
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$med_compensation,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Compensation",xlab="Dispersal rate",ylim=range(c(0,pmax(0,results$med_compensation,na.rm=TRUE))))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$med_compensation,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
  Corner_text(scenarios$network[i],"topright")
}
legend("right",scenarios$disturbance,col=dist_colours,lty=1,lwd=2,bty="n")
