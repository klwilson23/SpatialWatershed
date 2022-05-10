source("Functions/make networks.R")
source("Functions/Dispersal function.R")
source("Functions/patch_variance.R")
source("Functions/local disturbance.R")
source("Functions/some functions.R")
source("Functions/popDynFn.R")
source("Functions/finding MSY.R")
library(mvtnorm)
library(marima)
library(diagram)

metaPop <- function(Npatches=16,
                    networkType="linear",
                    patchDistance=1,
                    Nburnin=50,
                    NyrsPost=50,
                    omega=0.01,
                    m=100,
                    alpha=2,
                    metaK=1600,
                    cv=1e-6,
                    DistScenario="uniform",
                    magnitude_of_decline=0.9,
                    lagTime=1,
                    prodType="Beverton-Holt",
                    rho.time=1e-5,
                    rho.dist=1e5,
                    compensationLag=25,
                    dataWeighting=0.1,
                    alphaVariable=TRUE,
                    kVariable=TRUE,
                    spatialPlots=TRUE,
                    patchID=NA)
{
  log_sig <- sqrt(log(cv^2+1))
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
  
  MSY_p <- unlist(findMSY(alpha=alpha_p,beta=k_p,Npatches=Npatches,model=prodType)$MSY["yield",])
  RMSY_p <- unlist(findMSY(alpha=alpha_p,beta=k_p,Npatches=Npatches,model=prodType)$MSY["recruits",])
  NMSY_p <- unlist(findMSY(alpha=alpha_p,beta=k_p,Npatches=Npatches,model=prodType)$MSY["adults",])
  

  
  Nyears <- Nburnin+NyrsPost
  
  # part i - set some empty arrays
  
  # patch-specific dynamics to track
  popDyn <- array(NA,dim=c(Nyears,Npatches,2),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Stage"=c("Recruits","Spawners")))
  # ecological metrics to track over time:
  meta.rec <- rep(NA,Nyears)
  sink <- source <- pseudoSink <- var.rec <- array(NA,dim=c(Nyears,Npatches),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches))
  compensationBias <- rep(NA,Nyears)
  lostCapacity <- rep(NA,Nyears)
  alphaYr <- rep(NA,Nyears)
  metaKYr <- rep(NA,Nyears)
  MSYYr <- rep(NA,Nyears)
  surplusProd <- rep(NA,Nyears)
  Yield_MSY <- rep(NA,Nyears)
  NMSY_y <- rep(NA,Nyears)
  
  # dispersal information to track
  dispersing <- array(NA,dim=c(Nyears,Npatches,3),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Disersing"=c("Residents","Immigrants","Emigrants")))
  # metapopulation dynamics to track
  MetaPop <- matrix(NA,nrow=Nyears,ncol=3,dimnames=list("Year"=1:Nyears,"Stage"=c("Recruits","Spawners","Spawners_Contiguous")))
  # part ii - initialize populations
  rec.dev <- rnorm(Npatches,mean=0,sd=log_sig)
  meta.dev <- rnorm(1,mean=0,sd=log_sig)
  
  popDyn[1:lagTime,,"Spawners"] <- k_p
  popDyn[1:lagTime,,"Recruits"] <- k_p*exp(rec.dev-log_sig^2/2)
  MetaPop[1:lagTime,"Recruits"] <- sum(popDyn[1:lagTime,,"Recruits"])
  MetaPop[1:lagTime,"Spawners"] <- sum(popDyn[1:lagTime,,"Spawners"])
  MetaPop[1:lagTime,"Spawners_Contiguous"] <- metaK
  
  var.rec[1:lagTime,] <- rec.dev
  meta.rec[1:lagTime] <- meta.dev
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
    rec.dev <- dvSpaceTime(mnSig=log_sig^2,lastDV=var.rec[Iyear-lagTime,],rhoTime=rho.time,rhoSpa=rho.dist,distMatrix = distance_matrix) # mnSig is the variance terms for the variance-covariance matrix, adjusted on the scale of variance (not standard deviation)
    var.rec[Iyear,] <- rec.dev
    meta.rec[Iyear] <- rho.time*meta.rec[Iyear-1]+log_sig^2*(1-rho.time^2)
    rec.obs <- pmax(0,patch_rec*exp(rec.dev-log_sig^2/2))
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
    
    # part vi - calculate metapopulation dynamics & ecological metrics
    sink[Iyear,] <- (popDyn[Iyear,,"Recruits"] < popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
    source[Iyear,] <- (popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] > dispersing[Iyear,,"Immigrants"])
    pseudoSink[Iyear,] <- ((popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (popDyn[Iyear,,"Recruits"] > (popDyn[Iyear-lagTime,,"Spawners"]-dispersing[Iyear-lagTime,,"Immigrants"]))) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
    
    MetaPop[Iyear,"Recruits"] <- sum(popDyn[Iyear,,"Recruits"])
    MetaPop[Iyear,"Spawners"] <- sum(popDyn[Iyear,,"Spawners"])
    MetaPop[Iyear,"Spawners_Contiguous"] <- pmax(0,popDynamics(alpha=alpha,beta=metaK,Nadults=sum(popDyn[Iyear-lagTime,,"Spawners"]),model=prodType)$recruits*exp(meta.rec[Iyear]-log_sig^2/2))
    spawnRec <- data.frame("recruits"=MetaPop[(1+lagTime):Iyear,"Recruits"],"spawners"=MetaPop[1:(Iyear-lagTime),"Spawners"],weights=exp(dataWeighting*((1+lagTime):Iyear-Iyear))/max(exp(dataWeighting*((1+lagTime):Iyear-Iyear))))
    
    # part vii - resource assessment after disturbance
    if(Iyear > (Nburnin) & (MetaPop[Iyear,"Recruits"] > 0)) #Iyear > (Nburnin) & (MetaPop[Iyear,"Recruits"] > 0)
    {
      if(prodType=="Beverton-Holt"){
        alphaHat <- max(1.01,alphaLstYr,na.rm=TRUE) # initial values must be positive
        metaK_hat <- max(1,metaKLstYr,na.rm=TRUE)
        theta <- as.vector(c(alphaHat,log(metaK_hat),log(0.2)))

        SRfit <- try(optim(theta,SRfn,method="L-BFGS-B",lower=c(1.01,-Inf,-Inf),data=spawnRec,lastYr=list("alphaLstYr"=mean(c(alpha,alphaHat)),"metaKLstYr"=mean(c(metaK_hat,metaK)))), silent=T) # optimize likelihood function
        check<-is.numeric(SRfit[[1]]) # check to see if model converged
        
        ## store values only if model converged
        if(check[[1]] == "TRUE"){
          alphaHat <- SRfit$par[[1]]
          metaK_hat <- exp(SRfit$par[[2]])
          compHat <- (alphaHat*MetaPop[Iyear-1,"Spawners"])/(1+((alphaHat-1)/metaK_hat)*MetaPop[Iyear-1,"Spawners"])
          alphaLstYr <- alphaYr[Iyear] <- alphaHat
          metaKLstYr <- metaKYr[Iyear] <- metaK_hat
        }else{ # if model didn't converge, fit a Ricker
          SRfit <- lm(log(recruits/spawners)~spawners,data=spawnRec,weights=spawnRec$weights)
          alphaHat <- pmax(1.01,exp(coef(SRfit)["(Intercept)"]),na.rm=TRUE)
          metaK_hat <- pmax(1,-log(alphaHat)/coef(SRfit)[2],na.rm=TRUE)
          compHat <- alphaHat*MetaPop[Iyear-1,"Spawners"]*exp(-log(alphaHat)/metaK_hat*MetaPop[Iyear-1,"Spawners"])
          alphaLstYr <- alphaYr[Iyear] <- alphaHat
          metaKLstYr <- metaKYr[Iyear] <- metaK_hat
        }
      }else{ # if model truly is Ricker, then fit a Ricker
        SRfit <- lm(log(recruits/spawners)~spawners,data=spawnRec,weights=spawnRec$weights)
        alphaHat <- pmax(1.01,exp(coef(SRfit)["(Intercept)"]),na.rm=TRUE)
        metaK_hat <- pmax(1,-log(alphaHat)/coef(SRfit)[2],na.rm=TRUE)
        compHat <- alphaHat*MetaPop[Iyear-1,"Spawners"]*exp(-log(alphaHat)/metaK_hat*MetaPop[Iyear-1,"Spawners"])
        alphaLstYr <- alphaYr[Iyear] <- alphaHat
        metaKLstYr <- metaKYr[Iyear] <- metaK_hat
      } # estimate this year's estimated and realized production, MSY, etc.
      actualComp <- MetaPop[Iyear,"Recruits"]
      actualK <- sum(k_p)
      compensationBias[Iyear] <- compHat/actualComp
      lostCapacity[Iyear] <- metaK_hat/sum(k_p)
      
      tryMSY <- tryCatch(
        findMSY(alpha=alphaYr[Iyear],beta=metaKYr[Iyear],Npatches=1,model=prodType),
        error=function(e) e
      )
      
      if(!inherits(tryMSY, "error")){
        Yield_MSY[Iyear] <- unlist(tryMSY$MSY["yield",])
        NMSY_y[Iyear] <- unlist(tryMSY$MSY["adults",])
      }else{
        Yield_MSY[Iyear] <- 0
        NMSY_y[Iyear] <- 0
      }
      
      MSYYr[Iyear] <- Yield_MSY[Iyear]/sum(MSY_p)
      surplusProd[Iyear] <- (MetaPop[Iyear,"Spawners"]-NMSY_y[Iyear])/sum(popDyn[Iyear,,"Spawners"]-NMSY_p)
    }
    if(Iyear > (Nburnin) & (MetaPop[Iyear,"Recruits"] == 0))
    { # if the populaton went extinct, set estimated production and MSY to 0
      compHat <- 0
      metaK_hat <- 0
      actualComp <- MetaPop[Iyear,"Recruits"] # realized production
      actualK <- sum(k_p)
      compensationBias[Iyear] <- compHat/actualComp
      lostCapacity[Iyear] <- metaK_hat/sum(k_p)
      MSYYr[Iyear] <- 0/sum(MSY_p)
      surplusProd[Iyear] <- (MetaPop[Iyear,"Spawners"]-0)/sum(popDyn[Iyear,,"Spawners"]-NMSY_p)
    }
  }
  
  recovered <- half_recovered <- rep(FALSE,Nyears)
  recovered[(Nburnin+1):(Nyears-4)] <- running.mean(MetaPop[(Nburnin+1):Nyears,"Spawners"],5)>=mean(MetaPop[1:Nburnin,"Spawners"])
  
  half_recovered[(Nburnin+1):(Nyears-4)] <- running.mean(MetaPop[(Nburnin+1):Nyears,"Spawners"],5)>=(0.5*mean(MetaPop[1:Nburnin,"Spawners"]))
  
  recovery <- ifelse(any(recovered),min(which(recovered))-Nburnin,NyrsPost)
  half_recovery <- ifelse(any(half_recovered),min(which(half_recovered))-Nburnin,NyrsPost)
  
  extinction <- ifelse(any(MetaPop[,"Spawners"]==0),min(which(MetaPop[,"Spawners"]==0))-Nburnin,NyrsPost)
  patchOccupancy <- sapply(1:Nyears,function(x){sum(popDyn[x,,"Recruits"]>(0.1*k_p))/Npatches})
  
  distYear <- Nburnin+1 # what is the target reference year
  
  shortTermCV <- sd(popDyn[distYear:(Nburnin+half_recovery),,"Spawners"])/mean(popDyn[distYear:(Nburnin+half_recovery),,"Spawners"])
  medTermCV <- sd(popDyn[distYear:(Nburnin+recovery),,"Spawners"])/mean(popDyn[distYear:(Nburnin+recovery),,"Spawners"])
  longTermCV <- sd(popDyn[distYear:(Nyears),,"Spawners"])/mean(popDyn[distYear:(Nyears),,"Spawners"])
  
  shortTermProd <- mean(compensationBias[distYear:(distYear+5)],na.rm=TRUE)
  shortTermComp <- mean(alphaYr[distYear:(distYear+5)]/alpha,na.rm=TRUE)
  shortTermCap <- mean(metaKYr[distYear:(distYear+5)]/sum(k_p),na.rm=TRUE)
  shortTermMSY <- mean(MSYYr[distYear:(distYear+5)],na.rm=TRUE)
  shortTermOcc <- mean(patchOccupancy[distYear:(distYear+5)],na.rm=TRUE)
  shortTermSurp <- mean(MetaPop[distYear:(distYear+5),"Spawners"]/MetaPop[distYear:(distYear+5),"Spawners_Contiguous"],na.rm=TRUE)

  medTermProd <- mean(compensationBias[distYear:(distYear+10)],na.rm=TRUE)
  medTermComp <- mean(alphaYr[distYear:(distYear+10)]/alpha,na.rm=TRUE)
  medTermCap <- mean(metaKYr[distYear:(distYear+10)]/sum(k_p),na.rm=TRUE)
  medTermMSY <- mean(MSYYr[distYear:(distYear+10)],na.rm=TRUE)
  medTermOcc <- mean(patchOccupancy[distYear:(distYear+10)],na.rm=TRUE)
  medTermSurp <- mean(MetaPop[distYear:(distYear+10),"Spawners"]/MetaPop[distYear:(distYear+10),"Spawners_Contiguous"],na.rm=TRUE)
  
  longTermProd <- mean(compensationBias[distYear:(distYear+25)],na.rm=TRUE)
  longTermComp <- mean(alphaYr[distYear:(distYear+25)]/alpha,na.rm=TRUE)
  longTermCap <- mean(metaKYr[distYear:(distYear+25)]/sum(k_p),na.rm=TRUE)
  longTermMSY <- mean(MSYYr[distYear:(distYear+25)],na.rm=TRUE)
  longTermOcc <- mean(patchOccupancy[distYear:(distYear+25)],na.rm=TRUE)
  longTermSurp <- mean(MetaPop[distYear:(distYear+25),"Spawners"]/MetaPop[distYear:(distYear+25),"Spawners_Contiguous"],na.rm=TRUE)
  
  mnNMSY <- mean(NMSY_y[distYear:(Nburnin+recovery)],na.rm=TRUE)
  
  if(spatialPlots)
  {
    network
    spatialRecoveryPlot(textSize=1,popDyn,MetaPop,k_p,Nlevels=10,recovery,Nburnin,Nyears,alpha,metaK,alphaYr,metaKYr,lostCapacity,compensationBias,MSYYr,nodeScalar=35,network=network,networkType=networkType,Npatches=Npatches,NMsy=c(sum(NMSY_p),mnNMSY)) 
  }

  return(list("MetaPop"=MetaPop,"popDyn"=popDyn,"sink"=sink,"source"=source,"pseudoSink"=pseudoSink,"dispersing"=dispersing,"shortTermProd"=shortTermProd,"shortTermComp"=shortTermComp,"shortTermCap"=shortTermCap,"medTermProd"=medTermProd,"medTermComp"=medTermComp,"medTermCap"=medTermCap,"longTermProd"=longTermProd,"longTermComp"=longTermComp,"longTermCap"=longTermCap,"recovery"=recovery,"extinction"=extinction,"shortTermOcc"=shortTermOcc,"medTermOcc"=medTermOcc,"longTermOcc"=longTermOcc,"shortTermMSY"=shortTermMSY,"medTermMSY"=medTermMSY,"longTermMSY"=longTermMSY,"shortTermCV"=shortTermCV,"medTermCV"=medTermCV,"longTermCV"=longTermCV,"shortTermSurp"=shortTermSurp,"medTermSurp"=medTermSurp,"longTermSurp"=longTermSurp,"k_p"=k_p,"network"=network))
}
