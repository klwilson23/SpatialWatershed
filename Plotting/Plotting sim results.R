bias <- readRDS("bias.rds")
capacity <- readRDS("capacity.rds")
compensation <- readRDS("compensation.rds")
extinction <- readRDS("extinction.rds")
extinctionRate <- readRDS("extinctionRate.rds")
metaAbund <- readRDS("metaAbund.rds")
patchOcc <- readRDS("patchOcc.rds")
recovered <- readRDS("recovered.rds")
recovery <- readRDS("recovery.rds")


library(mvtnorm)
#library(marima)
require(mvtnorm)

source("Linear network.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("Metapop function.R")
source("popDynFn.R")

Nboots <- 50
Nlevels <- 5
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
alpha <- 3
metaK <- 1600
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.9
# what is the lag time between recruits and spawners
lagTime <- 1

network_levels <- c("linear","dendritic","star","complex")
disturbance_levels <- c("uniform","random","random_patch","targeted")
dispersal_levels <- c(0,seq(1e-3,0.1,length.out=Nlevels-1))
alpha_levels <- c("same","variable")
beta_levels <- c("same","variable")
spatial_levels <- c(100,10,1) # from low spatial dependency to high
temporal_levels <- c(1e-5,0.2,0.6) # from low temporal autocorrelation to high
stochastic_levels <- c(1e-1,0.2) # coefficient of variation on lognormal recruitment deviates


matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange","black")
for(i in 1:length(network_levels))
{
  matplot(dispersal_levels,t(recovery[i,,,"same","same",1,1,1]),col=dist_colours,lty=i,type="b",lwd=2,pch=21,ylab="Recovery",xlab="Dispersal rate",ylim=range(c(0,recovery[,,,"same","same",1,1,1])))
}

matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange","black")
for(i in 1:length(network_levels))
{
  matplot(dispersal_levels,t(recovery[i,,,"variable","same",1,1,1]),col=dist_colours,lty=i,type="b",lwd=2,pch=21,ylab="Recovery",xlab="Dispersal rate",ylim=range(c(0,recovery[,,,"variable","variable",1,1,1])))
}


matLayout <- matrix(1:4,nrow=2,ncol=2)
layout(matLayout)
dist_colours <- c("grey50","dodgerblue","orange","black")
for(i in 1:length(network_levels))
{
  matplot(dispersal_levels,t(capacity[i,,,"variable","variable",1,1,1]),col=dist_colours,lty=i,type="b",lwd=2,pch=21,ylab="Metapopulation capacity",xlab="Dispersal rate",ylim=range(c(0,capacity[,,,"variable","variable",1,1,1])))
}
