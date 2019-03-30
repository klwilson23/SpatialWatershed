
library(mvtnorm)
library(marima)

source("Linear network.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("Metapop function.R")

Nboots <- 100
Nlevels <- 10

# simulation parameters
Npatches <- 16
patchDist <- 1
Nburnin <- 50
Nyears <- Nburnin+100
m <- 1 # distance decay function: could set at some proportion of max distance
# adult stock-juvenile recruitment traits
alpha <- 3
metaK <- 1000
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.9
# what is the lag time between recruits and spawners
lagTime <- 1

network_levels <- c("linear","star","dendritic","complex")
disturbance_levels <- c("uniform","random","random_localized","targetted")
dispersal_levels <- c(0,exp(seq(log(1e-3),log(0.25),length.out=Nlevels-1)))
alpha_levels <- c("same","variable")
beta_levels <- c("same","variable")
spatial_levels <- exp(seq(log(1e-3),log(10),length.out=Nlevels))
temporal_levels <- seq(0,1,length.out=Nlevels)
stochastic_levels <- c(0.01,0.3,0.5,1.0)

recovery <- array(NA,dim=c(length(network_levels),
                           length(disturbance_levels),
                           length(dispersal_levels),
                           length(alpha_levels),
                           length(beta_levels),
                           length(spatial_levels),
                           length(temporal_levels),
                           length(stochastic_levels)),
                  dimnames=list("network"=network_levels,
                             "disturbance"=disturbance_levels,
                             "dispseral"=dispersal_levels,
                             "alpha"=alpha_levels,
                             "beta"=beta_levels,
                             "spatial"=spatial_levels,
                             "temporal"=temporal_levels,
                             "stochastic"=stochastic_levels))

length(recovery)

for(iNet in 1:length(network_levels))
  {
  makeNetworks(network_levels[iNet],Npatches=Npatches,patchDist=patchDist)
  for(iDisturb in 1:length(disturbance_levels))
    {
    for(iDisperse in 1:length(dispersal_levels))
      {
      for(iAlpha in 1:length(alpha_levels))
        {
        for(iBeta in 1:length(beta_levels))
          {
          for(iSpatial in 1:length(spatial_levels))
            {
            for(iTemporal in 1:length(temporal_levels))
              {
              for(iStochastic in 1:length(stochastic_levels))
                {
                recovery[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- SOMETHING
              }
            }
          }
        }
      }
    }
  }
}