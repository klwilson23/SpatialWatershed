
library(mvtnorm)
library(marima)

source("Linear network.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")


distance_matrix <- distances(linearLandscape,v=V(linearLandscape),to=V(linearLandscape))

ricker <-function(alpha,beta,Nadults){alpha*Nadults*exp(beta*Nadults)}

# leading parameters
# number of years & ecological scenarios
Nburnin <- 50
Nyears <- Nburnin+100
Npatches <- ncol(distance_matrix)
omega <- 0.01 # proportion of animals in patch that move
m <- 1 # distance decay function: could set at some proportion of max distance
# adult stock-juvenile recruitment traits
alpha <- 8
metaK <- 1000
beta <- -log(alpha)/metaK
cv <- 0.8
# temporal correlation
rho.time <- 0.7
# distance penalty to spatial correlation: higher means more independent
rho.dist <- 0.5
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.9
# what is the lag time between recruits and spawners
lagTime <- 1

alpha_heterogeneity <- TRUE
cap_heterogeneity <- TRUE
DistScenario <- "uniform"

# get the metapopulation & patch level Ricker parameters
alpha_p <- rep(alpha,Npatches)
k_p <- rep((metaK/Npatches),Npatches)
beta_p <- -log(alpha_p)/k_p
# call function to get the among patch variability in demographic traits
patches <- patch_variance(alpha_heterogeneity,cap_heterogeneity,Npatches,alpha_p,k_p)

alpha_p <- patches$alpha_p
beta_p <- patches$beta_p
k_p <- patches$k_p

# set some empty arrays

popDyn <- array(NA,dim=c(Nyears,Npatches,2),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Stage"=c("Recruits","Spawners")))

sink <- source <- psuedoSink <- var.rec <- array(NA,dim=c(Nyears,Npatches),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches))

dispersing <- array(NA,dim=c(Nyears,Npatches,3),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Disersing"=c("Residents","Immigrants","Emigrants")))
MetaPop <- matrix(NA,nrow=Nyears,ncol=2,dimnames=list("Year"=1:Nyears,"Stage"=c("Recruits","Spawners")))

# initialize populations
example.rec <- 1000*exp(rnorm(1e6,mean=0,sd=sqrt(log(cv^2+1)))-(log(cv^2+1))/2)
median(example.rec)
mean(example.rec)

rec.dev <- rnorm(Npatches,mean=0,sd=sqrt(log(cv^2+1)))

popDyn[1:lagTime,,"Spawners"] <- k_p
popDyn[1:lagTime,,"Recruits"] <- k_p*exp(rec.dev-(log(cv^2+1))/2)

MetaPop[1:lagTime,"Spawners"] <- sum(popDyn[1,,"Spawners"])

var.rec[1:lagTime,] <- rec.dev

for(Iyear in (lagTime+1):Nyears)
{
  
  # part iii - add disturbance
  if(Iyear==(Nburnin+1))
  {
    deaths_p <- Disturbance(metaPop=MetaPop[Iyear-1,"Spawners"],magnitude=magnitude_of_decline,DisType=DistScenario, N_p=popDyn[Iyear-1,,"Spawners"],prod=k_p)$deaths_p
    popDyn[Iyear-1,,"Spawners"] <- popDyn[Iyear-1,,"Spawners"]-deaths_p
    MetaPop[Iyear-1,"Spawners"] <- sum(popDyn[Iyear-1,,"Spawners"])
  }else{
    deaths_p <- rep(0,Npatches)
  }
  
  # part iia - population dynamics
  patch_rec <- ricker(alpha=alpha_p,beta=beta_p,popDyn[Iyear-lagTime,,"Spawners"])

  # part iib - stochastic recruitment
  
  rec.dev <- dvSpaceTime(mnSig=sqrt(log(cv^2+1)),lastDV=var.rec[Iyear-lagTime,],rhoTime=rho.time,rhoSpa=rho.dist,distMatrix = distance_matrix)
  var.rec[Iyear,] <- rec.dev
  
  rec.obs <- pmax(0,patch_rec*exp(rec.dev-(log(cv^2+1))/2))
  popDyn[Iyear,,"Recruits"] <- round(rec.obs)
  
  # part ii - dispersal between patches
  
  disperse <- dispersal(maxDispersal=omega,distDecay=m,dist_matrix=distance_matrix,recruitment=round(rec.obs))
  
  im <- disperse$immigrants
  em <- disperse$emigrants
  
  dispersing[Iyear,,"Residents"] <- popDyn[Iyear,,"Recruits"]-em
  dispersing[Iyear,,"Immigrants"] <- im
  dispersing[Iyear,,"Emigrants"] <- em
  
  popDyn[Iyear,,"Spawners"] <- popDyn[Iyear,,"Recruits"] + dispersing[Iyear,,"Immigrants"] - dispersing[Iyear,,"Emigrants"]
  
  # not sure these calculations should have recruitment lag times or not
  
  sink[Iyear,] <- (popDyn[Iyear,,"Recruits"] < popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
  
  source[Iyear,] <- (popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] > dispersing[Iyear,,"Immigrants"])
  
  psuedoSink[Iyear,] <- ((popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (popDyn[Iyear,,"Recruits"] > (popDyn[Iyear-lagTime,,"Spawners"]-dispersing[Iyear-lagTime,,"Immigrants"]))) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
  

  # part iv: calculate metapopulation dynamics
  MetaPop[Iyear,"Recruits"] <- sum(popDyn[Iyear,,"Recruits"])
  MetaPop[Iyear,"Spawners"] <- sum(popDyn[Iyear,,"Spawners"])
  
}

matplot(popDyn[,,"Spawners"]/k_p,type="l",xlab="Time",ylab="Relative abundance (N/K)")
lines(MetaPop[,"Spawners"]/metaK,lwd=3,col="black")
#return(list("patchDyn"=popDyn,"metaDyn"=MetaPop,"disDyn"=dispersing))
#plot(popDyn[1:(Nyears-1),20,"Spawners"],popDyn[2:Nyears,20,"Recruits"])
#plot(MetaPop[1:(Nyears-1),"Spawners"],MetaPop[2:Nyears,"Recruits"],type="p",xlab="Metapopulation spawners",ylab="Metapopulation recruits")

#plot(log(MetaPop[2:Nyears,"Recruits"]/MetaPop[1:(Nyears-1),"Spawners"]))

#matplot(dispersing[,,"Emigrants"],type="l")
#matplot(var.rec[,1:5],type="l")
#plot(var.rec[,1],type="l")
#acf(var.rec[,4])
MetaPop[50,"Spawners"]/MetaPop[49,"Spawners"]

