
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

source("Linear network.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("Metapop function.R")
source("popDynFn.R")
source("Functions/finding MSY.R")

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

results <- readRDS("Simulations/results2019-05-19.rds")
scenarios <- readRDS("Simulations/scenarios2019-05-19.rds")

NaNsims <- which(apply(apply(results[,-(1:8)],2,is.nan),1,any))

nScenarios <- length(NaNsims)

counter <- 0

progBar <- txtProgressBar(min = 0,  max = nScenarios, style = 3, initial = 0)
ptm <- Sys.time()

for(Isim in NaNsims)
{
  # call function to get the among patch variability in demographic traits
  alphaVariable <- ifelse(results$alpha[Isim]=="same",FALSE,TRUE)
  kVariable <- ifelse(results$beta[Isim]=="same",FALSE,TRUE)
  boatyMcboot <- replicate(Nboots,metaPop(Npatches=Npatches,
                                          networkType=results$network[Isim],
                                          patchDistance=patchDist,
                                          Nburnin=Nburnin,
                                          NyrsPost=NyrsPost,
                                          omega=results$dispersal[Isim],
                                          m=m,
                                          alpha=alpha,
                                          metaK=metaK,
                                          cv=results$stochastic[Isim],
                                          DistScenario=results$disturbance[Isim],
                                          magnitude_of_decline=magnitude_of_decline,
                                          lagTime=lagTime,
                                          prodType=model,
                                          rho.time=results$temporal[Isim],
                                          rho.dist=results$spatial[Isim],
                                          compensationLag=compLag,
                                          dataWeighting=dataWeighting,
                                          alphaVariable=alphaVariable,
                                          kVariable=kVariable,
                                          spatialPlots = FALSE))
  
  results$recovery[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="recovery"),]))
  
  results$recovered[Isim] <- sum(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="recovery"),])!=(NyrsPost))/Nboots
  
  results$extinctionRate[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="extinction"),]))
  
  results$extinction[Isim] <- sum(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="extinction"),])!=(NyrsPost))/Nboots
  
  results$shortOcc[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermOcc"),]))
  
  results$medOcc[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermOcc"),]))
  
  results$longOcc[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermOcc"),]))
  
  results$shortCV[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermCV"),]))
  
  results$medCV[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermCV"),]))
  
  results$longCV[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermCV"),]))
  
  results$shortMSY[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermMSY"),]))
  
  results$medMSY[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermMSY"),]))
  
  results$longMSY[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermMSY"),]))
  
  results$short_bias[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermProd"),]))
  
  results$short_capacity[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermCap"),]))
  
  results$short_compensation[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="shortTermComp"),]))
  results$med_bias[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermProd"),]))
  
  results$med_capacity[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermCap"),]))
  
  results$med_compensation[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="medTermComp"),]))
  
  results$long_bias[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermProd"),]))
  
  results$long_capacity[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermCap"),]))
  
  results$long_compensation[Isim] <- mean(unlist(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="longTermComp"),]))
  
  results$metaAbund[Isim] <- mean(unlist(lapply(boatyMcboot[which(dimnames(boatyMcboot)[[1]]=="MetaPop"),],function(x){x[Nburnin+NyrsPost,"Spawners"]})))/metaK
  counter <- counter + 1
  setTxtProgressBar(progBar, counter)
}
endtime <- Sys.time()-ptm
endtime

saveRDS(results,"Simulations/results2019-05-19_fixed.rds")
