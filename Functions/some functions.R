library(ggplot2)
Corner_text <- function(text, location="topright",...) #function to write text to the corner of plots
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}

wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")

get_beta <- function(mean,cv) #function that returns the alpha and beta shape parameters of a beta distribution, based on the mean and variation of a given beta distribution
{
  sd <- mean*cv
  alpha <- -((mean*(mean^2+sd^2-mean))/sd^2)
  beta <- alpha/mean-alpha
  return(list(alpha=alpha,beta=beta))
}

layout_in_circles <- function(g, group=1) {
  layout <- lapply(split(V(g), group), function(x) {
    layout_in_circle(induced_subgraph(g,x))
  })
  layout <- Map(`*`, layout, seq_along(layout))
  x <- matrix(0, nrow=vcount(g), ncol=2)
  split(x, group) <- layout
  x
}

extend <- function(alphabet) function(i) {
  base10toA <- function(n, A) {
    stopifnot(n >= 0L)
    N <- length(A)
    j <- n %/% N 
    if (j == 0L) A[n + 1L] else paste0(Recall(j - 1L, A), A[n %% N + 1L])
  }   
  vapply(i-1L, base10toA, character(1L), alphabet)
}

MORELETTERS <- extend(LETTERS)

dvSpaceTime <- function(mnSig,lastDV,rhoTime,rhoSpace,distMatrix)
{
  # function to calculate recruitment deviations correlated through space and time
  ar1 <- rhoTime*lastDV
  dvST <- ar1 + rmvnorm(1,mean=rep(0,ncol(distMatrix)),sigma=(mnSig*(1-rhoTime^2)*(exp(-rhoSpace*distMatrix))))
  return(dvST)
}

SRfn <- function(theta,data,lastYr){
  a.hat <- theta[1]
  b.hat <- exp(theta[2])
  sd.hat <- exp(theta[3])
  rec.mean <- (a.hat*data$spawners)/(1+((a.hat-1)/b.hat)*data$spawners)
  nll <- -1*sum(dlnorm(data$recruits,meanlog=log(rec.mean),sdlog=sd.hat,log=TRUE)*data$weights,na.rm=TRUE)
  penalty1 <- -dnorm(a.hat,lastYr$alphaLstYr,4*lastYr$alphaLstYr,log=TRUE) # penalized likelihood on estimated alpha
  penalty2 <- -dnorm(b.hat,lastYr$metaKLstYr,4*lastYr$metaKLstYr,log=TRUE) # penalized likelihood on estimated carrying capacity
  jnll <- sum(c(nll,penalty1,penalty2),na.rm=TRUE)
  return(jnll)
}

plottingFunc <- function(network,type,nodeScalar,Npatches,plotID=NA){
  par(xpd=TRUE,usr=c(-1.572561,1.572561,-1.080000,1.080000),mar=c(1,1,1,1))
  if(type=="linear"){
    plot(network$landscape,col="dodgerblue",layout=cbind(seq(-1,1,length.out = gorder(network$landscape)),seq(1,-1,length.out = gorder(network$landscape))),vertex.size=network$node.size*nodeScalar,xlim=c(-1,1),ylim=c(-1,1),rescale=FALSE,vertex.label=plotID)
  }
  if(type=="dendritic"){
    plot(network$landscape,col="dodgerblue",layout=layout_as_tree(network$landscape,root=V(network$landscape)[1]),vertex.size=network$node.size*nodeScalar,vertex.label=plotID)
  }
  if(type=="star"){
    plot(network$landscape,col="dodgerblue",layout=layout.reingold.tilford(network$landscape,circular=T),vertex.size=network$node.size*nodeScalar,vertex.label=plotID)
  }
  if(type=="complex"){
    
    spatialLayout <- matrix(MORELETTERS(1:Npatches),nrow=sqrt(Npatches),ncol=sqrt(Npatches),byrow=TRUE)
    
    spatialLayout <- t(sapply(1:Npatches,function(x){which(spatialLayout==LETTERS[x],arr.ind=TRUE)}))
    spatialLayout <- spatialLayout[rank(attr(V(network$landscape),"names")),]
    plot(network$landscape,col="dodgerblue",layout=spatialLayout,vertex.size=network$node.size*nodeScalar,vertex.label=plotID)
  }
}

spatialRecoveryPlot <- function(textSize=1,popDyn,MetaPop,k_p,Nlevels=10,recovery,Nburnin,Nyears,alpha,metaK,alphaYr,metaKYr,lostCapacity,compensationBias,MSY,nodeScalar=35,network,networkType=networkType,Npatches=Npatches,NMsy,patchID=NA)
{
  colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))
  
  matLayout <- matrix(0,ncol=18,nrow=18,byrow=T)
  matLayout[1:9,1:9] <- 1
  matLayout[1:4,10:12] <- 2
  matLayout[1:4,13:15] <- 3
  matLayout[5:8,10:12] <- 4
  matLayout[5:8,13:15] <- 5
  matLayout[9:12,10:12] <- 6
  matLayout[9:12,13:15] <- 7
  matLayout[4:9,16:17] <- 8
  matLayout[10:18,1:9] <- 9
  matLayout[13:18,10:18] <- 10
  layout(matLayout)
  par(mar=c(5,5,1,1))
  levelFactors <- factor(pmax(0.0,pmin(1.0,round(popDyn[Nyears,,"Spawners"]/k_p*Nlevels)/Nlevels)),levels=((0:10)/Nlevels))
  matplot(t(t(popDyn[,,"Spawners"])/k_p),type="l",xlab="Time",ylab="Relative abundance (N/K)",col=rev(colfunc(Nlevels+1))[levelFactors],ylim=c(0,1.1*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE)),cex.lab=textSize,xpd=TRUE)
  lines(MetaPop[,"Spawners"]/metaK,lwd=3,col="black")
  segments(y0=0,y1=1.05*max(t(popDyn[,,"Spawners"])/k_p),x0=(recovery+Nburnin),lwd=4,col="orange")
  
  recoverText <- ifelse((recovery+Nburnin)==Nyears,"Not recovered","Time to recovery")
  text(x=(recovery+Nburnin)-30,y=1.1*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE),recoverText,cex=textSize*0.8)
  curvedarrow(from=c((recovery+Nburnin)-30,1.08*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE)),to=c((recovery+Nburnin),1.0*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE)),lwd=2,lty=1,lcol="grey50",arr.col="grey50",curve=0.002,endhead=TRUE,arr.pos=0.65)
  
  text(x=(Nburnin)-25,y=0.05,"Time at ~90% loss",cex=textSize*0.8)
  curvedarrow(from=c(Nburnin-30,0.1),to=c(Nburnin+3,0.15),lwd=2,lty=1,lcol="grey50",arr.col="grey50",curve=-0.002,endhead=TRUE,arr.pos=0.65)
  
  for(i in c(Nburnin,Nburnin+3,Nburnin+6,Nburnin+9,Nburnin+12,Nburnin+15))
  {
    levelFactors <- factor(pmax(0.0,pmin(1.0,round(popDyn[i,,"Spawners"]/k_p*Nlevels)/Nlevels)),levels=((0:10)/Nlevels))
    V(network$landscape)$color <- rev(colfunc(Nlevels+1))[levelFactors]
    plottingFunc(network,networkType,nodeScalar,Npatches)
    title(main=paste("Year",i),line=0,font.main=1,cex.main=textSize)
  }

  par(mar=c(1,1,1,1))
  plot(NA,xlim=c(-1,1),ylim=c(-1,1),xaxt="n",yaxt="n",xlab="",ylab="",frame.plot=FALSE,xpd=NA)
  colorlegend(rev(colfunc(Nlevels)),zlevels=Nlevels,zlim=c(0,1),dz=0.1,posx=c(0.35,0.5),posy=c(0.1,0.95),digit=2)
  title(main=expression('N'[t]*'/K'),line=0.75,font.main=1,cex.main=textSize,xpd=NA)
  #text(x=0.5,y=1.75,expression('N'[t]*'/K'),cex=1.5,xpd=NA)
  
  par(mar=c(5,5,1,1))
  #plot(MetaPop[1:(Nyears-1),"Spawners"],MetaPop[2:Nyears,"Recruits"],type="p",xlab="Metapopulation adults",ylab="Metapopulation recruits",pch=21,bg=ifelse(1:(Nyears-1)>Nburnin,"orange","dodgerblue"),xlim=c(0,max(MetaPop[,"Spawners"],na.rm=TRUE)),ylim=c(0,max(MetaPop[,"Recruits"],na.rm=TRUE)),cex.lab=textSize)
  
  curve((alpha*x)/(1+((alpha-1)/metaK)*x),from=0,to=max(MetaPop[,"Spawners"],na.rm=TRUE),lwd=2,col="dodgerblue",xpd=FALSE,xlab="Metapopulation adults",ylab="Metapopulation recruits",cex.lab=textSize,ylim=c(0,max(MetaPop[,"Recruits"],na.rm=TRUE)))
  curve((mean(alphaYr[(Nburnin+1):(Nburnin+recovery)],na.rm=TRUE)*x)/(1+((mean(alphaYr[(Nburnin+1):(Nburnin+recovery)],na.rm=TRUE)-1)/mean(metaKYr[(Nburnin+1):(Nburnin+recovery)],na.rm=TRUE))*x),from=0,to=max(MetaPop[,"Spawners"],na.rm=TRUE),add=TRUE,lwd=2,col="orange",xpd=FALSE)
  
  #abline(v=NMsy[1],lwd=2,lty=2,col="dodgerblue",xpd=FALSE)
  #abline(v=NMsy[2],lwd=2,lty=2,col="orange",xpd=FALSE)
  
  
  legend("bottomright",c("Recrutment - pristine","Recruitment - disturbed"),pch=c(NA,NA),pt.bg=c("dodgerblue","orange"),lwd=c(1,1),lty=c(1,1),col=c("dodgerblue","orange"),bty="n",cex=textSize)
  
  par(mar=c(5,4,1,1))
  plot(lostCapacity[Nburnin:Nyears],xlab="Years after disturbance",ylab="Relative bias",type="l",ylim=range(c(-0.3,lostCapacity[Nburnin:Nyears],alphaYr[Nburnin:Nyears]/alpha,compensationBias[Nburnin:Nyears],MSY[Nburnin:Nyears]),na.rm=TRUE),lwd=2,col="dodgerblue",xpd=NA,cex.lab=textSize)
  lines(alphaYr[Nburnin:Nyears]/alpha,lwd=2,col="orange")
  lines(compensationBias[Nburnin:Nyears],lwd=2,col="grey50")
  lines(MSY[Nburnin:Nyears],lwd=2,col="black")
  legend("bottomright",c("capacity","compensation","production","MSY"),lty=1,col=c("dodgerblue","orange","grey50","black"),lwd=2,bty="n",cex=textSize)
}


spatialRecoveryPlotv2 <- function(textSize=1,popDyn,MetaPop,k_p,Nlevels=10,recovery,Nburnin,Nyears,nodeScalar=35,network,networkType=networkType,Npatches=Npatches)
{
  colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))
  par(mar=c(5,5,1,1))
  levelFactors <- factor(pmax(0.0,pmin(1.0,round(popDyn[Nyears,,"Spawners"]/k_p*Nlevels)/Nlevels)),levels=((0:10)/Nlevels))
  matplot(t(t(popDyn[,,"Spawners"])/k_p),type="l",xlab="Time",ylab="Relative abundance (N/K)",col=rev(colfunc(Nlevels+1))[levelFactors],ylim=c(0,1.1*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE)),cex.lab=textSize,xpd=TRUE)
  lines(MetaPop[,"Spawners"]/metaK,lwd=3,col="black")
  segments(y0=0,y1=1.05*max(t(popDyn[,,"Spawners"])/k_p),x0=(recovery+Nburnin),lwd=4,col="orange")
  
  recoverText <- ifelse((recovery+Nburnin)==Nyears,"Not recovered","Time at recovery")
  text(x=(recovery+Nburnin)-30,y=1.1*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE),recoverText,cex=textSize*0.8)
  curvedarrow(from=c((recovery+Nburnin)-30,1.08*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE)),to=c((recovery+Nburnin),1.0*max(t(popDyn[,,"Spawners"])/k_p,na.rm=TRUE)),lwd=2,lty=1,lcol="grey50",arr.col="grey50",curve=0.002,endhead=TRUE,arr.pos=0.65)
  
  text(x=(Nburnin)-25,y=0.05,"Time at ~90% loss",cex=textSize*0.8)
  curvedarrow(from=c(Nburnin-30,0.1),to=c(Nburnin+3,0.15),lwd=2,lty=1,lcol="grey50",arr.col="grey50",curve=-0.002,endhead=TRUE,arr.pos=0.65)
  
  for(i in c(Nburnin,Nburnin+3,Nburnin+6,Nburnin+9,Nburnin+12,Nburnin+15))
  {
    levelFactors <- factor(pmax(0.0,pmin(1.0,round(popDyn[i,,"Spawners"]/k_p*Nlevels)/Nlevels)),levels=((0:10)/Nlevels))
    V(network$landscape)$color <- rev(colfunc(Nlevels+1))[levelFactors]
    plottingFunc(network,networkType,nodeScalar,Npatches)
    title(main=paste("Year",i),line=0,font.main=1,cex.main=textSize)
  }
}

v_all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

findHull <- function(df,x,y)
{
  df[chull(df[,x],df[,y]),c(x,y)]
}

clustering <- function(x,n,type,standard)
{
  d <- dist(x,method="euclidean")
  
  layout(matrix(1,nrow=1,ncol=1))
  if(type=="hier")
  {
    fit <- hclust(d,method="ward.D2")
    boot.meth <- hclustCBI
    plot(fit)
  }
  if(type=="agnes")
  {
    fit <- agnes(d)
    boot.meth <- hclustCBI
    plot(fit)
  }
  if(type=="kmeans")
  {
    fit <- kmeans(d,n)
    boot.meth <- kmeansCBI
    
  }
  plot(fit)
  groups <- cutree(fit, k=n)
  rect.hclust(fit, k=n, border="red")
  clustStats <- cluster.stats(d,groups)
  SS <- clustStats$within.cluster.ss
  dunnFit <- clustStats$dunn # the Dunn index
  dunn <- clustStats$dunn2
  gamma <- clustStats$pearsongamma
  width <- clustStats$avg.silwidth
  widGap <- clustStats$widestgap
  separation <- clustStats$sindex

  return(list(data=d,fit=fit,groups=groups,dunnInd = dunnFit,
              dunn2nd = dunn,gamma=gamma,width=width,minSS=SS,gap=widGap,
              separat = separation))
}