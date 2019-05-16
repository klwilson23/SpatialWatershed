
# april 3, 2019:
# error happens for "targeted" disturbance where the simulation gets stuck in the while loop for accumulating patch disturbance
# solution? - if the sum of all high productivity patches is NOT enough to reach the TARGET LOSS then (1) take the next most productive patch or (2) nuke every last animal in that patch
# april 16, 2019:
# run the assessment for all available years of data after disturbance but with some data-weighting for recent years

# solutions: first 10 years post-disturbance is of interest (IUCN & COSEWIC criteria)

# april 23, 2019: we could use some level of "transience" or "variability" in the post-disturbance regime between year of disturbance and year of recovery.

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

Nboots <- 100
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
alpha <- 1.8
metaK <- 1600
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.9
# what is the lag time between recruits and spawners
lagTime <- 1

network_levels <- c("linear","dendritic","star","complex")
disturbance_levels <- c("uniform","random","random_patch")
dispersal_levels <- c(0,10^(seq(log10(1e-4),log10(0.1),length.out=Nlevels-1)))
alpha_levels <- c("same","variable")
beta_levels <- c("same","variable")
spatial_levels <- c(1e4,5,1e-1) # from low spatial dependency to high
temporal_levels <- c(1e-5,0.2,0.6) # from low temporal autocorrelation to high
stochastic_levels <- c(1e-4,1e-1) # coefficient of variation on lognormal recruitment deviates

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
                             "dispersal"=dispersal_levels,
                             "alpha"=alpha_levels,
                             "beta"=beta_levels,
                             "spatial"=spatial_levels,
                             "temporal"=temporal_levels,
                             "stochastic"=stochastic_levels))

shortCV <- medCV <- longCV <- shortMSY <- medMSY <- longMSY <- shortOcc <- medOcc <- longOcc <- metaAbund <- med_compensation <- med_capacity <- med_bias <- long_compensation <- long_capacity <- long_bias <- short_compensation <- short_capacity <- short_bias <- extinctionRate <- extinction <- recovered <- recovery

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
                # call function to get the among patch variability in demographic traits
                alphaVariable <- ifelse(alpha_levels[iAlpha]=="same",FALSE,TRUE)
                kVariable <- ifelse(beta_levels[iBeta]=="same",FALSE,TRUE)
                boatyMcboot <- replicate(Nboots,metaPop(Npatches=Npatches,
                                                        networkType=network_levels[iNet],
                                                        patchDistance=patchDist,
                                                        Nburnin=Nburnin,
                                                        NyrsPost=NyrsPost,
                                                        omega=dispersal_levels[iDisperse],
                                                        m=m,
                                                        alpha=alpha,
                                                        metaK=metaK,
                                                        cv=stochastic_levels[iStochastic],
                                                        DistScenario=disturbance_levels[iDisturb],
                                                        magnitude_of_decline=magnitude_of_decline,
                                                        lagTime=lagTime,
                                                        prodType=model,
                                                        rho.time=temporal_levels[iTemporal],
                                                        rho.dist=spatial_levels[iSpatial],
                                                        compensationLag=compLag,
                                                        dataWeighting=dataWeighting,
                                                        alphaVariable=FALSE,
                                                        kVariable=FALSE))
                
                recovery[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="recovery"),]))
                
                recovered[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- sum(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="recovery"),])!=(NyrsPost))
                
                extinctionRate[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="extinction"),]))
                
                extinction[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- sum(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="extinction"),])!=(NyrsPost))
                
                shortOcc[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermOcc"),]))
                
                medOcc[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermOcc"),]))
                
                longOcc[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermOcc"),]))
                
                shortCV[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermCV"),]))
                
                medCV[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermCV"),]))
                
                longCV[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermCV"),]))
                
                shortMSY[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermMSY"),]))
                
                medMSY[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermMSY"),]))
                
                longMSY[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermMSY"),]))
                
                short_bias[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermProd"),]))
                
                short_capacity[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermCap"),]))
                
                short_compensation[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermComp"),]))
                
                med_bias[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermProd"),]))
                
                med_capacity[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermCap"),]))
                
                med_compensation[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermComp"),]))
                
                long_bias[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermProd"),]))
                
                long_capacity[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermCap"),]))
                
                long_compensation[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermComp"),]))
                
                metaAbund[iNet,iDisturb,iDisperse,iAlpha,iBeta,iSpatial,iTemporal,iStochastic] <- mean(unlist(lapply(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="MetaPop"),],function(x){x[Nburnin+NyrsPost,"Spawners"]})))/metaK
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

saveRDS(recovery,paste("Simulations/recovery",Sys.Date(),".rds",sep=""))
saveRDS(recovered,paste("Simulations/recovered",Sys.Date(),".rds",sep=""))
saveRDS(extinction,paste("Simulations/extinction",Sys.Date(),".rds",sep=""))
saveRDS(extinctionRate,paste("Simulations/extinctionRate",Sys.Date(),".rds",sep=""))
saveRDS(shortOcc,paste("Simulations/shortOcc",Sys.Date(),".rds",sep=""))
saveRDS(medOcc,paste("Simulations/medOcc",Sys.Date(),".rds",sep=""))
saveRDS(longOcc,paste("Simulations/longOcc",Sys.Date(),".rds",sep=""))

saveRDS(shortCV,paste("Simulations/shortCV",Sys.Date(),".rds",sep=""))
saveRDS(medCV,paste("Simulations/medCV",Sys.Date(),".rds",sep=""))
saveRDS(longCV,paste("Simulations/longCV",Sys.Date(),".rds",sep=""))

saveRDS(shortMSY,paste("Simulations/shortMSY",Sys.Date(),".rds",sep=""))
saveRDS(medMSY,paste("Simulations/medMSY",Sys.Date(),".rds",sep=""))
saveRDS(longMSY,paste("Simulations/longMSY",Sys.Date(),".rds",sep=""))


saveRDS(short_bias,paste("Simulations/short_bias",Sys.Date(),".rds",sep=""))
saveRDS(short_compensation,paste("Simulations/short_compensation",Sys.Date(),".rds",sep=""))
saveRDS(short_capacity,paste("Simulations/short_capacity",Sys.Date(),".rds",sep=""))
saveRDS(med_bias,paste("Simulations/med_bias",Sys.Date(),".rds",sep=""))
saveRDS(med_compensation,paste("Simulations/med_compensation",Sys.Date(),".rds",sep=""))
saveRDS(med_capacity,paste("Simulations/med_capacity",Sys.Date(),".rds",sep=""))
saveRDS(long_bias,paste("Simulations/long_bias",Sys.Date(),".rds",sep=""))
saveRDS(long_compensation,paste("Simulations/long_compensation",Sys.Date(),".rds",sep=""))
saveRDS(long_capacity,paste("Simulations/long_capacity",Sys.Date(),".rds",sep=""))
saveRDS(metaAbund,paste("Simulations/metaAbund",Sys.Date(),".rds",sep=""))
