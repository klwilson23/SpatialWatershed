source("make networks.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("popDynFn.R")
library(mvtnorm)
library(marima)
library(diagram)

metaPop <- function(Npatches=16,
                    networkType="linear",
                    patchDistance=1,
                    Nburnin=50,
                    NyrsPost=100,
                    omega=0.01,
                    m=1,
                    alpha=2,
                    metaK=1600,
                    cv=1e-2,
                    DistScenario="uniform",
                    magnitude_of_decline=0.9,
                    lagTime=1,
                    prodType="Beverton-Holt",
                    rho.time=1e-5,
                    rho.dist=1e5,
                    compensationLag=25,
                    dataWeighting=0.1,
                    alphaVariable=FALSE,
                    kVariable=FALSE)
{
  # make the spatial network
  network <- makeNetworks(network=networkType,Npatches=Npatches,patchDist=patchDistance)
  distance_matrix <- network$distanceMatrix
  # get the metapopulation & patch level Ricker parameters
  alpha_p <- rep(alpha,Npatches)
  k_p <- rep((metaK/Npatches),Npatches)
  beta_p <- k_p
  # call function to get the among patch variability in demographic traits
  patches <- patch_variance(model=prodType,alphaVariable,kVariable,Npatches=Npatches,alpha=alpha,alpha_p=alpha_p,metaK=metaK,k_p=k_p,magnitude_of_decline=magnitude_of_decline)
  alpha_p <- patches$alpha_p
  beta_p <- patches$beta_p
  k_p <- patches$k_p
  
  Nyears <- Nburnin+NyrsPost
  
  # part i - set some empty arrays
  
  # patch-specific dynamics to track
  popDyn <- array(NA,dim=c(Nyears,Npatches,2),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Stage"=c("Recruits","Spawners")))
  # ecological metrics to track over time:
  sink <- source <- pseudoSink <- var.rec <- array(NA,dim=c(Nyears,Npatches),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches))
  compensationBias <- rep(NA,Nyears)
  lostCapacity <- rep(NA,Nyears)
  alphaYr <- rep(NA,Nyears)
  metaKYr <- rep(NA,Nyears)
  
  # dispersal information to track
  dispersing <- array(NA,dim=c(Nyears,Npatches,3),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Disersing"=c("Residents","Immigrants","Emigrants")))
  # metapopulation dynamics to track
  MetaPop <- matrix(NA,nrow=Nyears,ncol=2,dimnames=list("Year"=1:Nyears,"Stage"=c("Recruits","Spawners")))
  # part ii - initialize populations
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
      deaths_p <- Disturbance(metaPop=ceiling(MetaPop[Iyear-1,"Spawners"]),magnitude=magnitude_of_decline,DisType=DistScenario, N_p=ceiling(popDyn[Iyear-1,,"Spawners"]),prod=k_p,Npatches=Npatches)$deaths_p
      popDyn[Iyear-1,,"Spawners"] <- pmax(0,popDyn[Iyear-1,,"Spawners"]-deaths_p)
      MetaPop[Iyear-1,"Spawners"] <- sum(popDyn[Iyear-1,,"Spawners"])
    }else{
      deaths_p <- rep(0,Npatches)
    }
    
    # part iva - population dynamics
    patch_rec <- popDynamics(alpha=alpha_p,beta=beta_p,Nadults=popDyn[Iyear-lagTime,,"Spawners"],model=prodType)$recruits
    
    # part ivb - stochastic recruitment
    # function dvSpaceTime adds spatial & temporal correlation
    rec.dev <- dvSpaceTime(mnSig=sqrt(log(cv^2+1)),lastDV=var.rec[Iyear-lagTime,],rhoTime=rho.time,rhoSpa=rho.dist,distMatrix = distance_matrix)
    var.rec[Iyear,] <- rec.dev
    rec.obs <- pmax(0,patch_rec*exp(rec.dev-(log(cv^2+1))/2))
    popDyn[Iyear,,"Recruits"] <- rec.obs
    
    # part v - dispersal between patches
    # dispersal() function calculates emigration & immigration based on distance matrix and penalty to movement
    disperse <- dispersal(maxDispersal=omega,distDecay=m,dist_matrix=distance_matrix,recruitment=rec.obs)
    im <- disperse$immigrants
    em <- disperse$emigrants
    dispersing[Iyear,,"Residents"] <- popDyn[Iyear,,"Recruits"]-em
    dispersing[Iyear,,"Immigrants"] <- im
    dispersing[Iyear,,"Emigrants"] <- em
    popDyn[Iyear,,"Spawners"] <- pmax(0,popDyn[Iyear,,"Recruits"] + dispersing[Iyear,,"Immigrants"] - dispersing[Iyear,,"Emigrants"])
    
    # not sure these calculations should have recruitment lag times or not
    # part vi - calculate metapopulation dynamics & ecological metrics
    sink[Iyear,] <- (popDyn[Iyear,,"Recruits"] < popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
    source[Iyear,] <- (popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] > dispersing[Iyear,,"Immigrants"])
    pseudoSink[Iyear,] <- ((popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (popDyn[Iyear,,"Recruits"] > (popDyn[Iyear-lagTime,,"Spawners"]-dispersing[Iyear-lagTime,,"Immigrants"]))) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
    
    MetaPop[Iyear,"Recruits"] <- sum(popDyn[Iyear,,"Recruits"])
    MetaPop[Iyear,"Spawners"] <- sum(popDyn[Iyear,,"Spawners"])
    
    # part vii - estimate compensation
    if(Iyear > (Nburnin) & (MetaPop[Iyear,"Recruits"] > 0))
    {
      spawnRec <- data.frame("recruits"=MetaPop[2:Iyear,"Recruits"],"spawners"=MetaPop[1:(Iyear-1),"Spawners"],weights=exp(dataWeighting*(2:Iyear-Iyear))/max(exp(dataWeighting*(2:Iyear-Iyear))))
      
      if(prodType=="Beverton-Holt"){
        alphaHat <- pmax(1.01,alphaLstYr)
        metaK_hat <- pmax(1,metaKLstYr)
        theta <- as.vector(c(alphaHat,log(metaK_hat),log(0.2)))
        SRfit <- optim(theta,SRfn,method="L-BFGS-B",lower=c(1.01,-Inf,-Inf),data=spawnRec,lastYr=list("alphaLstYr"=alphaLstYr,"metaKLstYr"=metaKLstYr))
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
  extinction <- ifelse(any(MetaPop[,"Spawners"]==0),min(which(MetaPop[,"Spawners"]==0))-Nburnin,NyrsPost)
  patchOccupancy <- sum(popDyn[Nyears,,"Recruits"]>(0.05*k_p))/Npatches
  
  distYear <- Nburnin+1 # what is the target reference year
  shortTermProd <- mean(compensationBias[distYear:(distYear+5)],na.rm=TRUE)
  shortTermComp <- mean(alphaYr[distYear:(distYear+5)]/alpha,na.rm=TRUE)
  shortTermCap <- mean(metaKYr[distYear:(distYear+5)]/actualK,na.rm=TRUE)
  
  medTermProd <- mean(compensationBias[distYear:(distYear+10)],na.rm=TRUE)
  medTermComp <- mean(alphaYr[distYear:(distYear+10)]/alpha,na.rm=TRUE)
  medTermCap <- mean(metaKYr[distYear:(distYear+10)]/actualK,na.rm=TRUE)
  
  longTermProd <- mean(compensationBias[distYear:(distYear+25)],na.rm=TRUE)
  longTermComp <- mean(alphaYr[distYear:(distYear+25)]/alpha,na.rm=TRUE)
  longTermCap <- mean(metaKYr[distYear:(distYear+25)]/actualK,na.rm=TRUE)
  
  #spatialRecoveryPlot(textSize=1,popDyn,MetaPop,k_p,Nlevels=10,recovery,Nburnin,Nyears,alpha,metaK,alphaYr,metaKYr,lostCapacity,compensationBias,nodeScalar=35,network=network,networkType=networkType,Npatches=Npatches)
  
  return(list("MetaPop"=MetaPop,"popDyn"=popDyn,"sink"=sink,"source"=source,"pseudoSink"=pseudoSink,"dispersing"=dispersing,"shortTermProd"=shortTermProd,"shortTermComp"=shortTermComp,"shortTermCap"=shortTermCap,"medTermProd"=medTermProd,"medTermComp"=medTermComp,"medTermCap"=medTermCap,"longTermProd"=longTermProd,"longTermComp"=longTermComp,"longTermCap"=longTermCap,"recovery"=recovery,"extinction"=extinction,"patchOccupancy"=patchOccupancy))
}

#makeMeta <- metaPop(networkType = "linear", m=100)
#makeMeta <- metaPop(networkType = "dendritic", m=100)
#makeMeta <- metaPop(networkType = "star", m=100)
#makeMeta <- metaPop(networkType = "complex", m=100)
