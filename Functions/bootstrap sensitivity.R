
# april 3, 2019:
# error happens for "targeted" disturbance where the simulation gets stuck in the while loop for accumulating patch disturbance
# solution? - if the sum of all high productivity patches is NOT enough to reach the TARGET LOSS then (1) take the next most productive patch or (2) remove every last animal in that patch
# april 16, 2019:
# run the assessment for all available years of data after disturbance but with some data-weighting for recent years

# solutions: first 10 years post-disturbance is of interest (IUCN & COSEWIC criteria)

# april 23, 2019: we could use some level of "transience" or "variability" in the post-disturbance regime between year of disturbance and year of recovery.

library(mvtnorm)
#library(marima)
require(mvtnorm)

source("Functions/Linear network.R")
source("Functions/Dispersal function.R")
source("Functions/patch_variance.R")
source("Functions/local disturbance.R")
source("Functions/some functions.R")
source("Functions/Metapop function.R")
source("Functions/popDynFn.R")
source("Functions/finding MSY.R")

Nboots <- 100
boot_scen <- c(10,100,500,1000)
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

network_levels <- c("linear","dendritic","star","complex")
disturbance_levels <- c("uniform","random","random_patch")
dispersal_levels <- c(0,10^(seq(log10(1e-4),log10(0.05),length.out=Nlevels-1)))
alpha_levels <- c("same","variable")
beta_levels <- c("same","variable")
spatial_levels <- c(1e4,5,1e-1) # from low spatial dependency to high
temporal_levels <- c(1e-5,0.2,0.6) # from low temporal autocorrelation to high
stochastic_levels <- c(1e-3,1e-1) # coefficient of variation on lognormal recruitment deviates

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

nScenarios <- length(recovery) # how many scenarios are we simulating

counter <- 0
results <- as.data.frame(expand.grid(dimnames(recovery)))
results$stochastic <- as.numeric(as.character(results$stochastic))
results$spatial <- as.numeric(as.character(results$spatial))
results$temporal <- as.numeric(as.character(results$temporal))
results$dispersal <- as.numeric(as.character(results$dispersal))
results <- data.frame(results,
                      "recovery"=rep(NA,nrow(results)),
                      "recovered"=rep(NA,nrow(results)),
                      "extinction"=rep(NA,nrow(results)),
                      "extinctionRate"=rep(NA,nrow(results)),
                      "short_bias"=rep(NA,nrow(results)),
                      "short_capacity"=rep(NA,nrow(results)),
                      "short_compensation"=rep(NA,nrow(results)),
                      "long_bias"=rep(NA,nrow(results)),
                      "long_capacity"=rep(NA,nrow(results)),
                      "long_compensation"=rep(NA,nrow(results)),
                      "med_bias"=rep(NA,nrow(results)),
                      "med_capacity"=rep(NA,nrow(results)),
                      "med_compensation"=rep(NA,nrow(results)),
                      "metaAbund"=rep(NA,nrow(results)),
                      "longOcc"=rep(NA,nrow(results)),
                      "medOcc"=rep(NA,nrow(results)),
                      "shortOcc"=rep(NA,nrow(results)),
                      "longMSY"=rep(NA,nrow(results)),
                      "medMSY"=rep(NA,nrow(results)),
                      "shortMSY"=rep(NA,nrow(results)),
                      "longCV"=rep(NA,nrow(results)),
                      "medCV"=rep(NA,nrow(results)),
                      "shortCV"=rep(NA,nrow(results)),
                      "short_surp"=rep(NA,nrow(results)),
                      "med_surp"=rep(NA,nrow(results)),
                      "long_surp"=rep(NA,nrow(results)))
tail(results[results$network=="dendritic" & results$stochastic==1e-1,1:8],30)
boot_results <- results[rep(9478,length(boot_scen)), ]
boot_results$Nboots <- boot_scen
progBar <- txtProgressBar(min = 0,  max = nrow(boot_results), style = 3, initial = 0)
ptm <- Sys.time()

for(bootSim in 1:length(boot_scen)) # 1:nrow(results) #8997:9000
{
  Nboots <- boot_scen[bootSim]
  # call function to get the among patch variability in demographic traits
  alphaVariable <- ifelse(boot_results$alpha[bootSim]=="same",FALSE,TRUE)
  kVariable <- ifelse(boot_results$beta[bootSim]=="same",FALSE,TRUE)
  boatyMcboot <- replicate(Nboots,metaPop(Npatches=Npatches,
                                          networkType=boot_results$network[bootSim],
                                          patchDistance=patchDist,
                                          Nburnin=Nburnin,
                                          NyrsPost=NyrsPost,
                                          omega=boot_results$dispersal[bootSim],
                                          m=m,
                                          alpha=alpha,
                                          metaK=metaK,
                                          cv=boot_results$stochastic[bootSim],
                                          DistScenario=boot_results$disturbance[bootSim],
                                          magnitude_of_decline=magnitude_of_decline,
                                          lagTime=lagTime,
                                          prodType=model,
                                          rho.time=boot_results$temporal[bootSim],
                                          rho.dist=boot_results$spatial[bootSim],
                                          compensationLag=compLag,
                                          dataWeighting=dataWeighting,
                                          alphaVariable=alphaVariable,
                                          kVariable=kVariable,
                                          spatialPlots = FALSE))
  
  boot_results$recovery[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="recovery"),]))
  
  boot_results$recovered[bootSim] <- sum(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="recovery"),])!=(NyrsPost))/Nboots
  
  boot_results$extinctionRate[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="extinction"),]))
  
  boot_results$extinction[bootSim] <- sum(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="extinction"),])!=(NyrsPost))/Nboots
  
  boot_results$shortOcc[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermOcc"),]))
  
  boot_results$medOcc[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermOcc"),]))
  
  boot_results$longOcc[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermOcc"),]))
  
  boot_results$shortCV[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermCV"),]))
  
  boot_results$medCV[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermCV"),]))
  
  boot_results$longCV[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermCV"),]))
  
  boot_results$shortMSY[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermMSY"),]))
  
  boot_results$medMSY[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermMSY"),]))
  
  boot_results$longMSY[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermMSY"),]))
  
  boot_results$short_bias[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermProd"),]))
  
  boot_results$short_capacity[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermCap"),]))
  
  boot_results$short_compensation[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermComp"),]))
  boot_results$med_bias[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermProd"),]))
  
  boot_results$med_capacity[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermCap"),]))
  
  boot_results$med_compensation[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermComp"),]))
  
  boot_results$long_bias[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermProd"),]))
  
  boot_results$long_capacity[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermCap"),]))
  
  boot_results$long_compensation[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermComp"),]))
  
  boot_results$short_surp[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermSurp"),]))
  
  boot_results$med_surp[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermSurp"),]))
  
  boot_results$long_surp[bootSim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermSurp"),]))
  
  boot_results$metaAbund[bootSim] <- mean(unlist(lapply(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="MetaPop"),],function(x){x[Nburnin+NyrsPost,"Spawners"]})))/metaK
  counter <- counter + 1
  setTxtProgressBar(progBar, counter)
}
endtime <- Sys.time()-ptm
endtime

saveRDS(boot_results,paste("Simulations/bootstrap_results",Sys.Date(),".rds",sep=""))
