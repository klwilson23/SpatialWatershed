
source("Linear network.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")

distance_matrix <- distances(linearLandscape,v=V(linearLandscape),to=V(linearLandscape))

ricker <-function(alpha,beta,Nadults){alpha*Nadults*exp(beta*Nadults)}

Npatches <- ncol(distance_matrix)
# leading parameters
omega <- 0.001 # proportion of animals in patch that move
m <- 100 # distance decay function: could set at some proportion of max distance
# adult stock-juvenile recruitment traits
alpha <- 1.8
metaK <- 25000
beta <- -log(alpha)/metaK
cv <- 0.5
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

stock <- k_p
rec_mean <- round(sapply(1:Npatches,function(x){alpha_p[x]*stock[x]*exp(beta_p[x]*stock[x])})) # starting values for the patches

metaPopSize <- sum(rec_mean)

popDyn <- array(NA,dim=c(Nyears,Npatches,2),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Stage"=c("Recruits","Spawners")))

sink <- source <- psuedoSink <- array(NA,dim=c(Nyears,Npatches),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches))

dispersing <- array(NA,dim=c(Nyears,Npatches,3),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Disersing"=c("Residents","Immigrants","Emigrants")))
MetaPop <- matrix(NA,nrow=Nyears,ncol=2,dimnames=list("Year"=1:Nyears,"Stage"=c("Recruits","Spawners")))

popDyn[1,,"Spawners"] <- k_p
popDyn[1,,"Recruits"] <- k_p

MetaPop[1,"Spawners"] <- sum(popDyn[1,,"Spawners"])

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
  # part iib - stochastic recruitment
  rec.obs <- pmax(0,rnorm(Npatches,mean=patch_rec,sd=patch_rec*cv))
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

modularity(popDyn[,,"Spawners"])
?modularity
cluster::clusplot(popDyn[,,"Spawners"])
library(cluster)
plot(hclust(dist(popDyn[Nyears,,"Spawners"])))
cutree(hclust(dist(popDyn[Nyears,,"Spawners"])),k=4)
