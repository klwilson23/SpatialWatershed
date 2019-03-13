
library(mvtnorm)
library(marima)

source("Linear network.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")

distance_matrix <- distances(linearLandscape,v=V(linearLandscape),to=V(linearLandscape))

ricker <-function(alpha,beta,Nadults){alpha*Nadults*exp(beta*Nadults)}

Npatches <- ncol(distance_matrix)
# leading parameters
omega <- 1e-1 # proportion of animals in patch that move
m <- 1 # distance decay function: could set at some proportion of max distance
# adult stock-juvenile recruitment traits
alpha <- 2
metaK <- 25000
beta <- -log(alpha)/metaK
cv <- 0.5
# temporal correlation
rho.time <- 0.5
# distance penalty to spatial correlation: higher means more independent
rho.dist <- 0.1

alpha_heterogeneity <- TRUE
cap_heterogeneity <- FALSE

# number of years & ecological scenarios
Nburnin <- 50
Nyears <- Nburnin+100

# get the metapopulation & patch level Ricker parameters
alpha_p <- rep(alpha,Npatches)
k_p <- rep((metaK/Npatches),Npatches)
beta_p <- -log(alpha_p)/k_p


patches <- patch_variance(alpha_heterogeneity,cap_heterogeneity,Npatches,alpha_p,k_p)

alpha_p <- patches$alpha_p
beta_p <- patches$beta_p
k_p <- patches$k_p

popDyn <- array(NA,dim=c(Nyears,Npatches,2),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Stage"=c("Recruits","Spawners")))

sink <- source <- psuedoSink <- var.rec <- array(NA,dim=c(Nyears,Npatches),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches))

dispersing <- array(NA,dim=c(Nyears,Npatches,3),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Disersing"=c("Residents","Immigrants","Emigrants")))
MetaPop <- matrix(NA,nrow=Nyears,ncol=2,dimnames=list("Year"=1:Nyears,"Stage"=c("Recruits","Spawners")))

rec.dev <- rnorm(Npatches,mean=0,sd=k_p*cv)

popDyn[1,,"Spawners"] <- k_p
popDyn[1,,"Recruits"] <- k_p + rec.dev

MetaPop[1,"Spawners"] <- sum(popDyn[1,,"Spawners"])

var.rec[1,] <- rec.dev

for(Iyear in 2:Nyears)
{
  
  # part iii - add disturbance
  if(Iyear==(Nburnin+1))
  {
    deaths_p <- Disturbance(metaPop=MetaPop[Iyear-1,"Spawners"],magnitude=0.9,DisType="random", N_p=popDyn[Iyear-1,,"Spawners"],prod=alpha_p)$deaths_p
  }else{
    deaths_p <- rep(0,Npatches)
  }
  popDyn[Iyear-1,,"Spawners"] <- popDyn[Iyear-1,,"Spawners"]-deaths_p
  MetaPop[Iyear-1,"Spawners"] <- sum(popDyn[Iyear-1,,"Spawners"])
  
  # part iia - population dynamics
  patch_rec <- ricker(alpha=alpha_p,beta=beta_p,popDyn[Iyear-1,,"Spawners"])
  #patch_rec_t_1 <- ricker(alpha=alpha_p,beta=beta_p,popDyn[Iyear-2,,"Spawners"])
  
  # part iib - stochastic recruitment
  
  var.rec.temp <- (patch_rec*cv)
  var.time <- rho.time*var.rec[Iyear-1,]# + rnorm(Npatches,mean=0,sd=1)
  rec.dev <- rmvnorm(1,var.time,sigma=(cv*sum(patch_rec))*(1-rho.time)*(exp(-rho.dist*distance_matrix)))
  var.rec[Iyear,] <- rec.dev
  
  rec.obs <- pmax(0,patch_rec+rec.dev)
  popDyn[Iyear,,"Recruits"] <- round(rec.obs)
  
  # part ii - dispersal between patches
  
  disperse <- dispersal(maxDispersal=omega,distDecay=m,dist_matrix=distance_matrix,recruitment=round(rec.obs))
  
  im <- round(disperse$immigrants)
  em <- round(disperse$emigrants)
  
  dispersing[Iyear,,"Residents"] <- popDyn[Iyear,,"Recruits"]-im
  dispersing[Iyear,,"Immigrants"] <- im
  dispersing[Iyear,,"Emigrants"] <- em
  
  popDyn[Iyear,,"Spawners"] <- popDyn[Iyear,,"Recruits"] + dispersing[Iyear,,"Immigrants"] - dispersing[Iyear,,"Emigrants"]
  
  sink[Iyear,] <- (popDyn[Iyear,,"Recruits"] < popDyn[Iyear-1,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
  
  source[Iyear,] <- (popDyn[Iyear,,"Recruits"] > popDyn[Iyear-1,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] > dispersing[Iyear,,"Immigrants"])
  
  psuedoSink[Iyear,] <- ((popDyn[Iyear,,"Recruits"] > popDyn[Iyear-1,,"Spawners"]) & (popDyn[Iyear,,"Recruits"] > (popDyn[Iyear-1,,"Spawners"]-dispersing[Iyear-1,,"Immigrants"]))) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
  

  # part iv: calculate metapopulation dynamics
  MetaPop[Iyear,"Recruits"] <- sum(popDyn[Iyear,,"Recruits"])
  MetaPop[Iyear,"Spawners"] <- sum(popDyn[Iyear,,"Spawners"])
  
}

plot(MetaPop[,"Spawners"]/metaK,type="l")
matplot(popDyn[,,"Spawners"]/k_p,type="l")
lines(MetaPop[,"Spawners"]/metaK,lwd=3,col="black")
#return(list("patchDyn"=popDyn,"metaDyn"=MetaPop,"disDyn"=dispersing))
acf(popDyn[,Npatches,"Spawners"])
cor(popDyn[,c(1,2),"Spawners"])

plot(popDyn[1:(Nyears-1),20,"Spawners"],popDyn[2:Nyears,20,"Recruits"])

modularity(popDyn[,,"Spawners"])
?modularity
cluster::clusplot(popDyn[,,"Spawners"])
library(cluster)
plot(hclust(dist(popDyn[Nyears,,"Spawners"])))
cutree(hclust(dist(popDyn[Nyears,,"Spawners"])),k=4)
