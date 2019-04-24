
library(mvtnorm)
library(marima)
library(diagram)

source("make networks.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("popDynFn.R")

networkType <- "complex"

network <- makeNetworks(networkType,16,1)

plottingFunc <- function(network,type){
  if(type=="linear"){
    plot(network$landscape,col="dodgerblue",layout=cbind(0,seq(1,-1,length.out = gorder(network$landscape))),vertex.size=network$node.size*nodeScalar,xlim=c(-1,1),ylim=c(-1,1),rescale=FALSE)
  }
  if(type=="dendritic"){
    plot(network$landscape,col="dodgerblue",layout=layout_as_tree(network$landscape,root=V(network$landscape)[1]),vertex.size=network$node.size*nodeScalar)
  }
  if(type=="star"){
    plot(network$landscape,col="dodgerblue",layout=layout.reingold.tilford(network$landscape,circular=T),vertex.size=network$node.size*nodeScalar)
  }
  if(type=="complex"){
    
    spatialLayout <- matrix(MORELETTERS(1:Npatches),nrow=sqrt(Npatches),ncol=sqrt(Npatches),byrow=TRUE)
    
    spatialLayout <- t(sapply(1:Npatches,function(x){which(spatialLayout==LETTERS[x],arr.ind=TRUE)}))
    spatialLayout <- spatialLayout[rank(attr(V(network$landscape),"names")),]
    plot(network$landscape,col="dodgerblue",layout=spatialLayout,vertex.size=network$node.size*nodeScalar)
  }
}

colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))
Nlevels <- 10 #  how many levels for plotting spatial outcomes
nodeScalar <- 40
distance_matrix <- network$distanceMatrix

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
m <- 1 # distance decay function: could set at some proportion of max distance
# adult stock-juvenile recruitment traits
alpha <- 2
metaK <- Npatches*100
beta <- -log(alpha)/metaK
cv <- 0.1
# temporal correlation
rho.time <- 1e-5
# distance penalty to spatial correlation: higher means more independent
rho.dist <- 1e6
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.90
# what is the lag time between recruits and spawners
lagTime <- 1

alpha_heterogeneity <- TRUE
cap_heterogeneity <- TRUE
DistScenario <- "targeted"

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
  #patch_rec <- ricker(alpha=alpha_p,beta=beta_p,popDyn[Iyear-lagTime,,"Spawners"])
  patch_rec <- popDynamics(alpha=alpha_p,beta=beta_p,Nadults=popDyn[Iyear-lagTime,,"Spawners"],model=model)$recruits
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
  if(Iyear>(Nburnin+1) & (MetaPop[Iyear,"Recruits"] > 0))
  {
    spawnRec <- data.frame("recruits"=MetaPop[2:Iyear,"Recruits"],"spawners"=MetaPop[1:(Iyear-1),"Spawners"],weights=(1-sqrt(abs(2:Iyear-Iyear))/max(sqrt(abs(2:Iyear-Iyear)))))
    
    if(model=="Beverton-Holt"){
      SRfitTry <- lm(log(recruits/spawners)~spawners,data=spawnRec,weights=spawnRec$weights)
      alphaHat <- pmax(1.01,exp(coef(SRfitTry)["(Intercept)"]))
      metaK_hat <- pmax(10,-log(alphaHat)/coef(SRfitTry)[2])
      theta <- as.vector(c(alphaHat,log(metaK_hat),log(summary(SRfitTry)$sigma)))
      SRfit <- optim(theta,SRfn,method="L-BFGS-B",lower=c(1.01,-Inf,-Inf))
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
extinction <- ifelse(any(MetaPop[,"Spawners"]==0),which.min(MetaPop[,"Spawners"]==0),NA)
patchOccupancy <- sum(popDyn[Nyears,,"Recruits"]>(0.05*k_p))/Npatches
postDistBias <- sum(compensationBias[!is.na(compensationBias)])
lostCompensation <- alphaHat/alpha

matLayout <- matrix(0,ncol=18,nrow=18,byrow=T)
matLayout[1:9,1:9] <- 1
matLayout[1:3,10:12] <- 2
matLayout[1:3,13:15] <- 3
matLayout[4:6,10:12] <- 4
matLayout[4:6,13:15] <- 5
matLayout[7:9,10:12] <- 6
matLayout[7:9,13:15] <- 7
matLayout[3:7,16:17] <- 8
matLayout[10:18,1:9] <- 9
matLayout[10:18,10:18] <- 10
layout(matLayout)
par(mar=c(5,4,1,1))
levelFactors <- factor(pmax(0.0,pmin(1.0,round(popDyn[Nyears,,"Spawners"]/k_p*Nlevels)/Nlevels)),levels=((0:10)/Nlevels))
matplot(t(t(popDyn[,,"Spawners"])/k_p),type="l",xlab="Time",ylab="Relative abundance (N/K)",col=rev(colfunc(Nlevels+1))[levelFactors],ylim=c(0,1.1*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE)))
lines(MetaPop[,"Spawners"]/metaK,lwd=3,col="black")
segments(y0=0,y1=1.05*max(t(popDyn[,,"Spawners"])/k_p),x0=(recovery+Nburnin),lwd=4,col="orange")

text(x=(recovery+Nburnin)-30,y=1.1*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE),"Time to recovery")
curvedarrow(from=c((recovery+Nburnin)-30,1.08*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE)),to=c((recovery+Nburnin),1.0*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE)),lwd=2,lty=1,lcol="grey50",arr.col="grey50",curve=0.002,endhead=TRUE,arr.pos=0.65)

par(mar=c(1,1,1,1))
for(i in c(Nburnin,Nburnin+3,Nburnin+6,Nburnin+9,Nburnin+12,Nburnin+15))
{
  levelFactors <- factor(pmax(0.0,pmin(1.0,round(popDyn[i,,"Spawners"]/k_p*Nlevels)/Nlevels)),levels=((0:10)/Nlevels))
  V(network$landscape)$color <- rev(colfunc(Nlevels+1))[levelFactors]
  plottingFunc(network,networkType)
  title(main=paste("Year",i),line=0,font.main=1,cex=0.8)
}

par(mar=c(1,1,1,1))
plot(NA,xlim=c(-1,1),ylim=c(-1,1),xaxt="n",yaxt="n",xlab="",ylab="",frame.plot=FALSE,xpd=NA)
colorlegend(rev(colfunc(Nlevels)),zlevels=Nlevels,zlim=c(0,1),dz=0.1,posx=c(0.35,0.5),posy=c(0.1,0.95),digit=2)
title(main=expression('N'[t]*'/K'),line=1,font.main=1,cex=0.8,xpd=NA)
#text(x=0.5,y=1.75,expression('N'[t]*'/K'),cex=1.5,xpd=NA)

par(mar=c(5,4,1,1))
plot(MetaPop[1:(Nyears-1),"Spawners"],MetaPop[2:Nyears,"Recruits"],type="p",xlab="Metapopulation spawners",ylab="Metapopulation recruits",pch=21,bg=ifelse(1:(Nyears-1)>Nburnin,"orange","dodgerblue"))
legend("bottomright",c("pre-disturbance","post-disturbance"),pch=21,pt.bg=c("dodgerblue","orange"),bty="n")

par(mar=c(5,4,1,1))
plot(lostCapacity[(Nburnin+1):Nyears],xlab="Years",ylab="Post-disturbance carrying capacity",type="l",ylim=c(0,max(lostCapacity[(Nburnin+1):Nyears],na.rm=TRUE)),lwd=3)
#plot(compensationBias,xlab="Years",ylab="% bias in production")
#cumsum(compensationBias[!is.na(compensationBias)])
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
