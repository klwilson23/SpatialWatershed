
library(mvtnorm)
library(marima)

source("Linear network.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("popDynFn.R")

distance_matrix <- distances(linearLandscape,v=V(linearLandscape),to=V(linearLandscape))

ricker <-function(alpha,beta,Nadults){alpha*Nadults*exp(beta*Nadults)}

# leading parameters
# number of years & ecological scenarios
Nburnin <- 50
Nyears <- Nburnin+100
compensationLag <- 25 # how many years to lag estimates of compensation
Npatches <- ncol(distance_matrix)
omega <- 1e-6 # proportion of animals in patch that move
m <- 1 # distance decay function: could set at some proportion of max distance
# adult stock-juvenile recruitment traits
alpha <- 3
metaK <- 1000
beta <- -log(alpha)/metaK
cv <- 0.1
# temporal correlation
rho.time <- 0.1
# distance penalty to spatial correlation: higher means more independent
rho.dist <- 100
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.9
# what is the lag time between recruits and spawners
lagTime <- 1

alpha_heterogeneity <- FALSE
cap_heterogeneity <- FALSE
DistScenario <- "random"

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
# patch-specific dynamics to track
popDyn <- array(NA,dim=c(Nyears,Npatches,2),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Stage"=c("Recruits","Spawners")))

# ecological metrics to track over time:
sink <- source <- psuedoSink <- var.rec <- array(NA,dim=c(Nyears,Npatches),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches))
compensationBias <- rep(NA,Nyears)
lostCapacity <- rep(NA,Nyears)

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
  # ecological metrics time
  
  sink[Iyear,] <- (popDyn[Iyear,,"Recruits"] < popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
  
  source[Iyear,] <- (popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] > dispersing[Iyear,,"Immigrants"])
  
  psuedoSink[Iyear,] <- ((popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (popDyn[Iyear,,"Recruits"] > (popDyn[Iyear-lagTime,,"Spawners"]-dispersing[Iyear-lagTime,,"Immigrants"]))) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
  

  # part iv: calculate metapopulation dynamics
  MetaPop[Iyear,"Recruits"] <- sum(popDyn[Iyear,,"Recruits"])
  MetaPop[Iyear,"Spawners"] <- sum(popDyn[Iyear,,"Spawners"])
  
  # part iii - add disturbance
  # part vii - estimate compensation
  if(Iyear>(compensationLag+1) & (MetaPop[Iyear,"Recruits"] > 0))
  {
    spawnRec <- data.frame("recruits"=MetaPop[((Iyear-compensationLag)+1):Iyear,"Recruits"],"spawners"=MetaPop[(Iyear-compensationLag):(Iyear-1),"Spawners"],weights=sqrt(((Iyear-compensationLag)+1):Iyear)/max(sqrt(((Iyear-compensationLag)+1):Iyear)))
    
    if(model=="Beverton-Holt"){
      SRfit <- nls(recruits~(a.hat*spawners)/(1+((a.hat-1)/b.hat)*spawners),data=spawnRec,weights=spawnRec$weights,start=list("a.hat"=3,"b.hat"=1000),lower=c(0,0),upper=c(30,Inf),algorithm="port")
      alphaHat <- coef(SRfit)["a.hat"]
      metaK_hat <- coef(SRfit)["b.hat"]
      compHat <- (alphaHat*MetaPop[Iyear-1,"Spawners"])/(1+((alphaHat-1)/metaK_hat)*MetaPop[Iyear-1,"Spawners"])
    }else{
      SRfit <- lm(log(recruits/spawners)~spawners,data=spawnRec,weights=spawnRec$weights)
      alphaHat <- exp(coef(SRfit)["(Intercept)"])
      metaK_hat <- -log(compHat)/coef(SRfit)[2]
      compHat <- alphaHat*MetaPop[Iyear-1,"Spawners"]*exp(-log(alphaHat)/metaK_hat*MetaPop[Iyear-1,"Spawners"])
    }
    actualComp <- MetaPop[Iyear,"Recruits"] # realized production
    actualK <- sum(k_p)
    compensationBias[Iyear] <- 100*(as.numeric(alphaHat-alpha)/alpha)
    lostCapacity[Iyear] <- 100*(as.numeric(metaK_hat-actualK)/actualK)
  }
}

layout(matrix(1:2,ncol=1,nrow=2,byrow=T))
par(mar=c(5,4,1,1))
matplot(popDyn[,,"Spawners"]/k_p,type="l",xlab="Time",ylab="Relative abundance (N/K)")
lines(MetaPop[,"Spawners"]/metaK,lwd=3,col="black")
#return(list("patchDyn"=popDyn,"metaDyn"=MetaPop,"disDyn"=dispersing))
#plot(popDyn[1:(Nyears-1),20,"Spawners"],popDyn[2:Nyears,20,"Recruits"])
plot(MetaPop[1:(Nyears-1),"Spawners"],MetaPop[2:Nyears,"Recruits"],type="p",xlab="Metapopulation spawners",ylab="Metapopulation recruits",pch=21,bg=ifelse(1:(Nyears-1)>Nburnin,"orange","dodgerblue"))
legend("topright",c("pre-disturbance","post-disturbance"),pch=21,pt.bg=c("dodgerblue","orange"),bty="n")

par(mar=c(5,4,1,1))
plot(lostCapacity,xlab="Years",ylab="% bias in carrying capacity")
plot(compensationBias,xlab="Years",ylab="% bias in production")

#plot(log(MetaPop[2:Nyears,"Recruits"]/MetaPop[1:(Nyears-1),"Spawners"]))

#matplot(dispersing[,,"Emigrants"],type="l")
#matplot(var.rec[,1:5],type="l")
#plot(var.rec[,1],type="l")
#acf(var.rec[,4])
MetaPop[50,"Spawners"]/MetaPop[49,"Spawners"]
round(MetaPop[49,"Spawners"]*(1-magnitude_of_decline))

# policy metrics
# years to recovery (averaged over 10 years)
# estimated compensation to true compensation 
# patch occupancy
