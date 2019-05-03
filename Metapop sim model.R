
library(mvtnorm)
library(marima)
library(diagram)

source("make networks.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("popDynFn.R")

networkType <- "linear"

network <- makeNetworks(networkType,16,1)

colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))
Nlevels <- 10 #  how many levels for plotting spatial outcomes
nodeScalar <- 40
distance_matrix <- network$distanceMatrix

dataWeighting <- 0.1#layout(1)
#curve(1+exp(dataWeighting*x)/max(exp(dataWeighting*x)),from=0,to=-100)

ricker <-function(alpha,beta,Nadults){alpha*Nadults*exp(beta*Nadults)}

# leading parameters
# number of years & ecological scenarios
model <- "Beverton-Holt"
Nburnin <- 50
NyrsPost <- 100
Nyears <- Nburnin+NyrsPost
targetYear <- Nburnin+10
compensationLag <- 50 # how many years to lag estimates of compensation
Npatches <- ncol(distance_matrix)
omega <- 0.01 # proportion of animals in patch that move
m <- 100 # distance decay function: could set at some proportion of max distance
# adult stock-juvenile recruitment traits
alpha <- 2
metaK <- Npatches*100
beta <- -log(alpha)/metaK
cv <- 1e-5
# temporal correlation
rho.time <- 1e-6
# distance penalty to spatial correlation: higher means more independent
rho.dist <- 1e6
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.90
# what is the lag time between recruits and spawners
lagTime <- 1

alpha_heterogeneity <- TRUE
cap_heterogeneity <- TRUE
DistScenario <- "random_patch"

# get the metapopulation & patch level Ricker parameters
alpha_p <- rep(alpha,Npatches)
k_p <- rep((metaK/Npatches),Npatches)
beta_p <- -log(alpha_p)/k_p
# call function to get the among patch variability in demographic traits
patches <- patch_variance(model=model,alpha_heterogeneity,cap_heterogeneity,Npatches,alpha=alpha,alpha_p=alpha_p,metaK=metaK,k_p=k_p,magnitude_of_decline=magnitude_of_decline)

alpha_p <- patches$alpha_p
beta_p <- patches$beta_p
k_p <- patches$k_p

# set some empty arrays
# patch-specific dynamics to track
popDyn <- array(NA,dim=c(Nyears,Npatches,2),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Stage"=c("Recruits","Spawners")))

# ecological metrics to track over time:
sink <- source <- psuedoSink <- var.rec <- array(NA,dim=c(Nyears,Npatches),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches))
compensationBias <- rep(NA,Nyears)
lostCapacity <- rep(NA,Nyears)
alphaYr <- rep(NA,Nyears)
metaKYr <- rep(NA,Nyears)


# dispersal information to track
dispersing <- array(NA,dim=c(Nyears,Npatches,3),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Disersing"=c("Residents","Immigrants","Emigrants")))

# metapopulation dynamics to track
MetaPop <- matrix(NA,nrow=Nyears,ncol=2,dimnames=list("Year"=1:Nyears,"Stage"=c("Recruits","Spawners")))

# initialize populations
rec.dev <- rnorm(Npatches,mean=0,sd=sqrt(log(cv^2+1)))

popDyn[1:lagTime,,"Spawners"] <- k_p
popDyn[1:lagTime,,"Recruits"] <- k_p*exp(rec.dev-(log(cv^2+1))/2)

MetaPop[1:lagTime,"Recruits"] <- sum(popDyn[1:lagTime,,"Recruits"])
MetaPop[1:lagTime,"Spawners"] <- sum(popDyn[1:lagTime,,"Spawners"])

var.rec[1:lagTime,] <- rec.dev

alphaLstYr <- alphaYr[1] <- alpha
metaKLstYr <- metaKYr[1] <- metaK

for(Iyear in (lagTime+1):Nyears)
{
  
  # part iii - add disturbance
  if(Iyear==(Nburnin+1))
  {
    deaths_p <- Disturbance(metaPop=ceiling(MetaPop[Iyear-1,"Spawners"]),magnitude=magnitude_of_decline,DisType=DistScenario, N_p=ceiling(popDyn[Iyear-1,,"Spawners"]),prod=k_p)$deaths_p
    popDyn[Iyear-1,,"Spawners"] <- pmax(0,popDyn[Iyear-1,,"Spawners"]-deaths_p)
    MetaPop[Iyear-1,"Spawners"] <- sum(popDyn[Iyear-1,,"Spawners"])
  }else{
    deaths_p <- rep(0,Npatches)
  }
  
  # part iia - population dynamics
  #patch_rec <- ricker(alpha=alpha_p,beta=beta_p,popDyn[Iyear-lagTime,,"Spawners"])
  patch_rec <- popDynamics(alpha=alpha_p,beta=beta_p,Nadults=popDyn[Iyear-lagTime,,"Spawners"],model=model)$recruits
  # part iib - stochastic recruitment
  
  rec.dev <- dvSpaceTime(mnSig=sqrt(log(cv^2+1)),lastDV=var.rec[Iyear-lagTime,],rhoTime=rho.time,rhoSpa=rho.dist,distMatrix = distance_matrix)
  var.rec[Iyear,] <- rec.dev
  
  rec.obs <- pmax(0,patch_rec*exp(rec.dev-(log(cv^2+1))/2))
  popDyn[Iyear,,"Recruits"] <- rec.obs
  
  # part ii - dispersal between patches
  
  disperse <- dispersal(maxDispersal=omega,distDecay=m,dist_matrix=distance_matrix,recruitment=rec.obs)
  
  im <- disperse$immigrants
  em <- disperse$emigrants
  
  dispersing[Iyear,,"Residents"] <- popDyn[Iyear,,"Recruits"]-em
  dispersing[Iyear,,"Immigrants"] <- im
  dispersing[Iyear,,"Emigrants"] <- em
  
  popDyn[Iyear,,"Spawners"] <- pmax(0,popDyn[Iyear,,"Recruits"] + dispersing[Iyear,,"Immigrants"] - dispersing[Iyear,,"Emigrants"])
  
  # not sure these calculations should have recruitment lag times or not
  # ecological metrics time
  
  sink[Iyear,] <- (popDyn[Iyear,,"Recruits"] < popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
  
  source[Iyear,] <- (popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] > dispersing[Iyear,,"Immigrants"])
  
  psuedoSink[Iyear,] <- ((popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (popDyn[Iyear,,"Recruits"] > (popDyn[Iyear-lagTime,,"Spawners"]-dispersing[Iyear-lagTime,,"Immigrants"]))) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
  

  # part iv: calculate metapopulation dynamics
  MetaPop[Iyear,"Recruits"] <- sum(popDyn[Iyear,,"Recruits"])
  MetaPop[Iyear,"Spawners"] <- sum(popDyn[Iyear,,"Spawners"])
  
  # part vii - estimate compensation
  if(Iyear > (Nburnin) & MetaPop[Iyear,"Recruits"] > 0)
  {
    spawnRec <- data.frame("recruits"=MetaPop[2:Iyear,"Recruits"],"spawners"=MetaPop[1:(Iyear-1),"Spawners"],weights=exp(dataWeighting*(2:Iyear-Iyear))/max(exp(dataWeighting*(2:Iyear-Iyear))))
    
    if(model=="Beverton-Holt"){
      SRfitTry <- lm(log(recruits/spawners)~spawners,data=spawnRec)
      alphaHat <- pmax(1.01,alphaLstYr)
      metaK_hat <- pmax(1,metaKLstYr)
      theta <- as.vector(c(alphaHat,log(metaK_hat),log(summary(SRfitTry)$sigma)))
      SRfit <- optim(theta,SRfn,method="L-BFGS-B",lower=c(1.01,-Inf,-Inf))
      alphaHat <- SRfit$par[1]
      metaK_hat <- exp(SRfit$par[2])
      compHat <- (alphaHat*MetaPop[Iyear-1,"Spawners"])/(1+((alphaHat-1)/metaK_hat)*MetaPop[Iyear-1,"Spawners"])
      alphaLstYr <- alphaYr[Iyear] <- alphaHat
      metaKLstYr <- metaKYr[Iyear] <- metaK_hat
    }else{
      SRfit <- lm(log(recruits/spawners)~spawners,data=spawnRec,weights=spawnRec$weights)
      alphaHat <- exp(coef(SRfit)["(Intercept)"])
      metaK_hat <- -log(alphaHat)/coef(SRfit)[2]
      compHat <- alphaHat*MetaPop[Iyear-1,"Spawners"]*exp(-log(alphaHat)/metaK_hat*MetaPop[Iyear-1,"Spawners"])
      alphaLstYr <- alphaYr[Iyear] <- alphaHat
      metaKLstYr <- metaKYr[Iyear] <- metaK_hat
    }
    actualComp <- MetaPop[Iyear,"Recruits"] # realized production
    actualK <- sum(k_p)
    compensationBias[Iyear] <- compHat/actualComp
    lostCapacity[Iyear] <- metaK_hat/actualK
  }
}

recovered <- rep(FALSE,Nyears)
recovered[(Nburnin+1):(Nyears-4)] <- running.mean(MetaPop[(Nburnin+1):Nyears,"Spawners"],5)>=mean(MetaPop[1:Nburnin,"Spawners"])
recovery <- ifelse(any(recovered),min(which(recovered))-Nburnin,NyrsPost)
extinction <- ifelse(any(MetaPop[,"Spawners"]==0),which.min(MetaPop[,"Spawners"]==0),NA)
patchOccupancy <- sum(popDyn[Nyears,,"Recruits"]>(0.05*k_p))/Npatches
postDistBias <- sum(compensationBias[!is.na(compensationBias)])
lostCompensation <- alphaHat/alpha

#jpeg(paste(networkType," example.jpeg",sep=""),res=800,units="in",height=8,width=9)
#nodeScalar=35
#spatialRecoveryPlot(textSize=0.8)
#dev.off()

cumsum(compensationBias[!is.na(compensationBias)]-1)
#plot(log(MetaPop[2:Nyears,"Recruits"]/MetaPop[1:(Nyears-1),"Spawners"]))

#matplot(dispersing[,,"Emigrants"],type="l")
#matplot(var.rec[,1:5],type="l")
#plot(var.rec[,1],type="l")
#acf(var.rec[,4])
MetaPop[50,"Spawners"]/MetaPop[49,"Spawners"]
round(MetaPop[49,"Spawners"]*(1-magnitude_of_decline))

comparisonYr <- 51
(alphaYr[comparisonYr]*MetaPop[comparisonYr-1,"Spawners"])/(1+((alphaYr[comparisonYr]-1)/metaKYr[comparisonYr])*MetaPop[comparisonYr-1,"Spawners"])
sum((alpha_p*popDyn[comparisonYr-1,,"Spawners"])/(1+((alpha_p-1)/k_p)*popDyn[comparisonYr-1,,"Spawners"]))
MetaPop[comparisonYr,"Recruits"]

# policy metrics
# years to recovery (averaged over 10 years)
# estimated compensation to true compensation 
# patch occupancy
