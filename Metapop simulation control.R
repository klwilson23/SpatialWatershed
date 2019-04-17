
# april 3, 2019:
# error happens for "targeted" disturbance where the simulation gets stuck in the while loop for accumulating patch disturbance
# solution? - if the sum of all high productivity patches is NOT enough to reach the TARGET LOSS then (1) take the next most productive patch or (2) nuke every last animal in that patch
# april 16, 2019:
# run the assessment for all available years of data after disturbance but with some data-weighting for recent years

# solutions: first 10 years post-disturbance is of interest (IUCN & COSEWIC criteria)

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
m <- 2 # distance decay function: penalty of 1 says about 60% of dispersing recruits move to neighbor patches. penalty of 2 says about 90% of dispersing recruits move to neighbor patches
# adult stock-juvenile recruitment traits
alpha <- 5
metaK <- 1600
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.9
# what is the lag time between recruits and spawners
lagTime <- 1

network_levels <- c("linear","star","dendritic","complex")
disturbance_levels <- c("uniform","random","random_patch","targeted")
dispersal_levels <- c(0,seq(1e-3,0.25,length.out=Nlevels-1))
alpha_levels <- c("same","variable")
beta_levels <- c("same","variable")
spatial_levels <- c(10,5,1) # from low spatial dependency to high
temporal_levels <- c(0.01,0.2,0.5) # from low temporal autocorrelation to high
stochastic_levels <- c(1e-6) # coefficient of variation on lognormal recruitment deviates

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

metaAbund <- compensation <- capacity <- bias <- patchOcc <- extinctionRate <- extinction <- recovered <- recovery

nScenarios <- length(recovery) # how many scenarios are we simulating

counter <- 0
#iNet = 1
#iDisturb=1
#iDisperse=1
#iAlpha=1
#iBeta=1
#iSpatial=1
#iTemporal=1
#iStochastic=1

progBar <- txtProgressBar(min = 0,  max = nScenarios, style = 3, initial = 0)
ptm = Sys.time()
for(iNet in 1:length(network_levels))
  {
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
                network <- makeNetworks(network_levels[iNet],Npatches=Npatches,patchDist=patchDist)
                # get the metapopulation & patch level Ricker parameters
                alpha_p <- rep(alpha,Npatches)
                k_p <- rep((metaK/Npatches),Npatches)
                beta_p <- k_p
                # call function to get the among patch variability in demographic traits
                alphaVariable <- ifelse(alpha_levels[iAlpha]=="same",FALSE,TRUE)
                kVariable <- ifelse(beta_levels[iBeta]=="same",FALSE,TRUE)
                patches <- patch_variance(model=model,alphaVariable,kVariable,Npatches,alpha=alpha,alpha_p=alpha_p,metaK=metaK,k_p=k_p)
                
                alpha_p <- patches$alpha_p
                beta_p <- patches$beta_p
                k_p <- patches$k_p
                if(iStochastic==1){ # if stochasticity is low, don't bootstrap
                  boatyMcboot <- replicate(1,metaPop(Npatches=Npatches,
                                                      distance_matrix=network$distanceMatrix,
                                                      Nburnin=Nburnin,
                                                      NyrsPost=NyrsPost,
                                                      omega=dispersal_levels[iDisperse],
                                                      m=m,
                                                      alpha=alpha,
                                                      metaK=metaK,
                                                      alpha_p=alpha_p,
                                                      beta_p=beta_p,
                                                      k_p=k_p,
                                                      cv=stochastic_levels[iStochastic],
                                                      DistScenario=disturbance_levels[iDisturb],
                                                      magnitude_of_decline=magnitude_of_decline,
                                                      lagTime=lagTime,
                                                      prodType=model,
                                                      rho.time=temporal_levels[iTemporal],
                                                      rho.dist=spatial_levels[iSpatial],
                                                      compensationLag=compLag))
                }else{
                  boatyMcboot <- replicate(Nboots,metaPop(Npatches=Npatches,
                                                          distance_matrix=network$distanceMatrix,
                                                          Nburnin=Nburnin,
                                                          NyrsPost=NyrsPost,
                                                          omega=dispersal_levels[iDisperse],
                                                          m=m,
                                                          alpha=alpha,
                                                          metaK=metaK,
                                                          alpha_p=alpha_p,
                                                          beta_p=beta_p,
                                                          k_p=k_p,
                                                          cv=stochastic_levels[iStochastic],
                                                          DistScenario=disturbance_levels[iDisturb],
                                                          magnitude_of_decline=magnitude_of_decline,
                                                          lagTime=lagTime,
                                                          prodType=model,
                                                          rho.time=temporal_levels[iTemporal],
                                                          rho.dist=spatial_levels[iSpatial],
                                                          compensationLag=compLag))
                }
                
                recovery[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="recovery"),]))
                
                recovered[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- sum(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="recovery"),])!=(NyrsPost))
                
                extinctionRate[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="extinction"),]))
                
                extinction[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- sum(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="extinction"),])!=(NyrsPost))
                
                patchOcc[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="patchOccupancy"),]))
                
                bias[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="bias"),]))
                
                capacity[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="lostCapacity"),]))
                
                compensation[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="lostCompensation"),]))
                
                metaAbund[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(lapply(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="MetaPop"),],function(x){x[Nburnin+NyrsPost,"Spawners"]})))/metaK
                #MetaPop <- boatyMcboot[[which(dimnames(boatyMcboot)[[1]]=="MetaPop")]]
                #popDyn <- boatyMcboot[[which(dimnames(boatyMcboot)[[1]]=="popDyn")]]
                #layout(matrix(1:2,ncol=1,nrow=2,byrow=T))
                #par(mar=c(5,4,1,1))
                #matplot(popDyn[,,"Spawners"]/k_p,type="l",xlab="Time",ylab="Relative abundance (N/K)")
                #lines(MetaPop[,"Spawners"]/metaK,lwd=3,col="black")
                #plot(MetaPop[1:((Nburnin+NyrsPost)-1),"Spawners"],MetaPop[2:(Nburnin+NyrsPost),"Recruits"],type="p",xlab="Metapopulation spawners",ylab="Metapopulation recruits",pch=21,bg=ifelse(1:((Nburnin+NyrsPost)-1)>Nburnin,"orange","dodgerblue"))
                #legend("topright",c("pre-disturbance","post-disturbance"),pch=21,pt.bg=c("dodgerblue","orange"),bty="n")
                #layout(1)
                counter <- counter + 1
                setTxtProgressBar(progBar, counter)
              }
            }
          }
        }
      }
    }
  }
}
endtime <- Sys.time()-ptm
endtime
