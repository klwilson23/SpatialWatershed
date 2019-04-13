
metaPop <- function(Npatches=16,distance_matrix,Nburnin=50,NyrsPost=100,omega=0.02,m=1,alpha=3,metaK=1000,alpha_p=rep(3,Npatches),beta_p=rep(1000,Npatches),k_p=rep(1000,Npatches),cv=0.01,DistScenario="uniform",magnitude_of_decline=0.9,lagTime=1,prodType="Beverton-Holt",rho.time=0,rho.dist=1e3,compensationLag=25)
{
  Nyears <- Nburnin+NyrsPost
  # part i - set some empty arrays
  # patch-specific dynamics to track
  popDyn <- array(NA,dim=c(Nyears,Npatches,2),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches,"Stage"=c("Recruits","Spawners")))
  # ecological metrics to track over time:
  sink <- source <- pseudoSink <- var.rec <- array(NA,dim=c(Nyears,Npatches),dimnames=list("Year"=1:Nyears,"Patch No."=1:Npatches))
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
    
    # part iva - population dynamics
    patch_rec <- popDynamics(alpha=alpha_p,beta=beta_p,Nadults=popDyn[Iyear-lagTime,,"Spawners"],model=model)$recruits
    
    # part ivb - stochastic recruitment
    # function dvSpaceTime adds spatial & temporal correlation
    rec.dev <- dvSpaceTime(mnSig=sqrt(log(cv^2+1)),lastDV=var.rec[Iyear-lagTime,],rhoTime=rho.time,rhoSpa=rho.dist,distMatrix = distance_matrix)
    var.rec[Iyear,] <- rec.dev
    rec.obs <- pmax(0,patch_rec*exp(rec.dev-(log(cv^2+1))/2))
    popDyn[Iyear,,"Recruits"] <- round(rec.obs)
    
    # part v - dispersal between patches
    # dispersal() function calculates emigration & immigration based on distance matrix and penalty to movement
    disperse <- dispersal(maxDispersal=omega,distDecay=m,dist_matrix=distance_matrix,recruitment=round(rec.obs))
    im <- disperse$immigrants
    em <- disperse$emigrants
    dispersing[Iyear,,"Residents"] <- popDyn[Iyear,,"Recruits"]-em
    dispersing[Iyear,,"Immigrants"] <- im
    dispersing[Iyear,,"Emigrants"] <- em
    popDyn[Iyear,,"Spawners"] <- popDyn[Iyear,,"Recruits"] + dispersing[Iyear,,"Immigrants"] - dispersing[Iyear,,"Emigrants"]
    
    # not sure these calculations should have recruitment lag times or not
    # part vi - calculate metapopulation dynamics & ecological metrics
    sink[Iyear,] <- (popDyn[Iyear,,"Recruits"] < popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
    source[Iyear,] <- (popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (dispersing[Iyear,,"Emigrants"] > dispersing[Iyear,,"Immigrants"])
    pseudoSink[Iyear,] <- ((popDyn[Iyear,,"Recruits"] > popDyn[Iyear-lagTime,,"Spawners"]) & (popDyn[Iyear,,"Recruits"] > (popDyn[Iyear-lagTime,,"Spawners"]-dispersing[Iyear-lagTime,,"Immigrants"]))) & (dispersing[Iyear,,"Emigrants"] < dispersing[Iyear,,"Immigrants"])
    
    MetaPop[Iyear,"Recruits"] <- sum(popDyn[Iyear,,"Recruits"])
    MetaPop[Iyear,"Spawners"] <- sum(popDyn[Iyear,,"Spawners"])
    
    # part vii - estimate compensation
    if(Iyear>(compensationLag+1) & (MetaPop[Iyear,"Recruits"] > 0))
    {
      spawnRec <- data.frame("recruits"=MetaPop[((Iyear-compensationLag)+1):Iyear,"Recruits"],"spawners"=MetaPop[(Iyear-compensationLag):(Iyear-1),"Spawners"],weights=sqrt(1:compensationLag)/max(sqrt(1:compensationLag)))
      if(model=="Beverton-Holt"){
        theta <- as.vector(c(alpha,log(metaK),log(1)))
        SRfit <- optim(theta,SRfn,method="L-BFGS-B",lower=c(1.01,0,-Inf),upper=c(50,10*metaK,Inf))
        alphaHat <- SRfit$par[1]
        metaK_hat <- exp(SRfit$par[2])
        compHat <- (alphaHat*MetaPop[Iyear-1,"Spawners"])/(1+((alphaHat-1)/metaK_hat)*MetaPop[Iyear-1,"Spawners"])
      }else{
        SRfit <- lm(log(recruits/spawners)~spawners,data=spawnRec,weights=spawnRec$weights)
        alphaHat <- exp(coef(SRfit)["(Intercept)"])
        metaK_hat <- -log(alphaHat)/coef(SRfit)[2]
        compHat <- alphaHat*MetaPop[Iyear-1,"Spawners"]*exp(-log(alphaHat)/metaK_hat*MetaPop[Iyear-1,"Spawners"])
      }
      actualComp <- MetaPop[Iyear,"Recruits"] # realized production
      actualK <- sum(k_p)
      compensationBias[Iyear] <- 100*(as.numeric(compHat-actualComp)/actualComp)
      lostCapacity[Iyear] <- metaK_hat/actualK
    }
  }
  recovered <- rep(FALSE,Nyears)
  recovered[(Nburnin+1):(Nyears-4)] <- running.mean(MetaPop[(Nburnin+1):Nyears,"Spawners"],5)>=mean(MetaPop[1:Nburnin,"Spawners"])
  recovery <- ifelse(any(recovered),min(which(recovered))-Nburnin,NyrsPost)
  extinction <- ifelse(any(MetaPop[,"Spawners"]==0),min(which(MetaPop[,"Spawners"]==0))-Nburnin,NyrsPost)
  patchOccupancy <- sum(popDyn[Nyears,,"Recruits"]>(0.05*k_p))/Npatches
  postDistBias <- sum(compensationBias[!is.na(compensationBias)])
  lostCompensation <- alphaHat/alpha
  
  return(list("MetaPop"=MetaPop,"popDyn"=popDyn,"sink"=sink,"source"=source,"pseudoSink"=pseudoSink,"dispersing"=dispersing,"bias"=postDistBias,"lostCompensation"=lostCompensation,"lostCapacity"=lostCapacity[Nyears],"recovery"=recovery,"extinction"=extinction,"patchOccupancy"=patchOccupancy))
}