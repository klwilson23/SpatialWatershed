results <- readRDS("Simulations/results2019-06-30.rds")
scenarios <- readRDS("Simulations/scenarios2019-06-30.rds")

library(mvtnorm)
library(vioplot)
library(ggplot2)
library(ggalluvial)
library("alphahull")
library(wesanderson)
library(ggpubr)

wesAnderson <- "Moonrise1"

source("Functions/Linear network.R")
source("Functions/Dispersal function.R")
source("Functions/patch_variance.R")
source("Functions/local disturbance.R")
source("Functions/some functions.R")
source("Functions/Metapop function.R")
source("Functions/popDynFn.R")
source("Functions/Add landscapes plot.R")

NaNsims <- which(apply(apply(results[,-(1:8)],2,is.nan),1,any))

results[NaNsims,]

Nboots <- 100
Nlevels <- 11
model <- "Beverton-Holt"
# simulation parameters
Npatches <- 16
patchDist <- 1
Nburnin <- 50
NyrsPost <- 100
Nyears <- Nburnin+NyrsPost
compLag <- 25 # how many years to lag estimates of compensation
dataWeighting <- 0.1 # penalty to past years for data weighting
m <- 100 # distance decay function: penalty of 1 says about 60% of dispersing recruits move to neighbor patches. penalty of 2 says about 90% of dispersing recruits move to neighbor patches
# adult stock-juvenile recruitment traits
alpha <- 2
metaK <- 1600
# how big is the disturbance after Nburnin years?
magnitude_of_decline <- 0.9
# what is the lag time between recruits and spawners
lagTime <- 1
results$StateShift <- ((1-results$recovered)+(1-results$longOcc))/2
results$RecoveryRate <- 1-results$recovery/NyrsPost
results$collapsed <- 1-results$recovered

results$unrecovered <- factor(ifelse(results$metaAbund >=0.9,"Recovered",ifelse(results$metaAbund>=0.7,"Recovering","Collapsed")),levels=c("Recovered","Recovering","Collapsed"))
results$recovery_pace <- factor(ifelse(results$recovery <=10,"Fast",ifelse(results$recovery<=25,"Moderate","Slow")),levels=c("Fast","Moderate","Slow"))
results$contraction <- factor(ifelse(results$medOcc >=0.9,"Full",ifelse(results$medOcc >=0.15,"Filling","Contracted")),levels=c("Full","Filling","Contracted"))
results$masked <- factor(ifelse(results$unrecovered!="Collapsed" & results$longOcc <= 0.50,"Hidden collapses","Good monitoring"),levels=c("Good monitoring","Hidden collapses"))
results$lostCap <- factor(ifelse(results$longMSY <=0.5,"Lost capacity","Healthy"),levels=c("Healthy","Lost capacity"))
results$outcome <- factor(paste(results$unrecovered,results$recovery_pace,results$contraction,results$masked,results$lostCap))
results$surprise <- factor(ifelse(grepl("Contracted",results$outcome),"Spatial contraction",
                           ifelse(grepl("Hidden",results$outcome),"Hidden collapses",
                                  ifelse(grepl("Lost capacity",results$outcome),"Lost capacity",
                                         ifelse(grepl("Fast",results$outcome) & grepl("Recover",results$outcome),"Resilient","Slow recovery")))),levels=c("Resilient","Slow recovery","Lost capacity","Hidden collapses","Spatial contraction"))

results$surprise_logic <- ifelse(results$surprise=="Resilient",1,1)
results$dispersal_range <- factor(ifelse(results$dispersal>=0.001,"High","Low"),levels=c("High","Low"))
results$network_lab <- factor(results$network,levels=levels(results$network),labels=c("Linear","Dendritic","Star","Complex"))
results$disturb_lab <- factor(results$disturbance,levels=levels(results$disturbance),labels=c("Uniform","Local, random","Local, extirpation"))
results$density_dep <- factor(results$alpha,levels=levels(results$alpha),labels=c("Identical","Diverse"))

surprises <- data.frame(aggregate(surprise_logic~network_lab+dispersal_range+disturb_lab+density_dep+surprise,data=results,FUN=sum))

gradColour <- colorRampPalette(rev(c("#bd0026","tomato","#fdae61","dodgerblue","forestgreen")))
margins <- c(0.1,1,0.1,0.1)
p1 <- ggplot(surprises,aes(y = surprise_logic,axis1 = density_dep, axis2 = dispersal_range,axis3=disturb_lab,axis4=surprise)) +
      geom_alluvium(aes(fill = surprise),width = 1/6,reverse=FALSE) +
      guides(fill = FALSE) +
      geom_stratum(width = 1/6,fill="grey",color="white", reverse = FALSE) +
      geom_text(stat = "stratum",size=3.5, infer.label = TRUE, reverse = FALSE) +
      scale_x_discrete(limits = c("Patch productivity", "Dispersal","Disturbance","Surprise"),expand=c(0.05,0.05)) +
      scale_fill_manual(values=gradColour(n=length(unique(surprises$surprise)))) +
      #scale_fill_brewer(type="div",palette="RdYlBu",direction=-1) +
      theme_minimal() +
      ylab("Frequency of outcome") +
      coord_cartesian(clip = "off") +
      theme(legend.position="none",strip.text.x=element_blank(),plot.margin=unit(margins,"line")) +
      ggtitle("Ecological surprises in metapopulations")
pAnnotated <- annotate_figure(p1,bottom=text_grob(wrapper("Interplay between density-dependent productivity, disturbance, and dispersal can lead to suprising recovery outcomes.",width=115),color="black",hjust=0,x=0.01,face="italic",size=10))
pAnnotated
wid_height <- 16/12
height <- 14
ggsave("Figures/surprising outcomes.jpeg",pAnnotated,units="cm",dpi=800,width=height*wid_height,height=height)

dist_colours <- c("tomato","dodgerblue","orange")
transparency <- 0.6
plotColours <- ifelse(results$dispersal==-1,
                           adjustcolor(dist_colours[results$disturbance],1,offset=c(-0.35,-0.35,-0.35,0)),
                           adjustcolor(dist_colours[results$disturbance],1))

tiff("Figures/Disturbance impacts on recovery regime polygon.tiff",compression="lzw",units="in",height=8,width=8,res=800)
layout(matrix(1:4,nrow=2,ncol=2,byrow=TRUE))
par(mar=c(4,4,0.5,0.5))
plot(results$RecoveryRate,results$longOcc,pch=21,bg=plotColours,xlab="Recovery rate (yr-1)",ylab="Patch occupancy (25 years)")

invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="RecoveryRate",y="longOcc"),col=adjustcolor(dist_colours[z],transparency))}))

plot(results$RecoveryRate,results$longMSY,pch=21,bg=plotColours,xlab="Recovery rate (yr-1)",ylab="Surplus production (25 years)")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="RecoveryRate",y="longMSY"),col=adjustcolor(dist_colours[z],transparency))}))

plot(results$RecoveryRate,results$collapsed,pch=21,bg=plotColours,xlab="Recovery rate (yr-1)",ylab="Collapse rate (% of simulations)")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="RecoveryRate",y="collapsed"),col=adjustcolor(dist_colours[z],transparency))}))

plot(results$longMSY,results$longOcc,pch=21,bg=plotColours,xlab="Surplus production (25 years)",ylab="Patch occupancy (25 years)")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="longMSY",y="longOcc"),col=adjustcolor(dist_colours[z],transparency))}))

legend("bottomright",c("Uniform","Local, random","Local, extirpation"),pch=22,pt.bg=dist_colours,bty="n",title="Disturbance regime")
dev.off()

tiff("Figures/Disturbance impacts on recovery regime smoothed polygon.tiff",compression="lzw",units="in",height=8,width=8,res=800)
layout(matrix(1:4,nrow=2,ncol=2,byrow=TRUE))
par(mar=c(4,4,0.5,0.5))
plot(results$RecoveryRate,results$longOcc,pch=21,bg=plotColours,xlab="Recovery rate (yr-1)",ylab="Patch occupancy (25 years)")

invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="RecoveryRate",y="longOcc"),col=adjustcolor(dist_colours[z],transparency))}))

plot(results$RecoveryRate,results$longMSY,pch=21,bg=plotColours,xlab="Recovery rate (yr-1)",ylab="Surplus production (25 years)")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="RecoveryRate",y="longMSY"),col=adjustcolor(dist_colours[z],transparency))}))

plot(results$RecoveryRate,results$collapsed,pch=21,bg=plotColours,xlab="Recovery rate (yr-1)",ylab="Collapse rate (% of simulations)")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){
  if(z<3){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="RecoveryRate",y="collapsed"),col=adjustcolor(dist_colours[z],transparency))}else{xMat <- unique(cbind(results[results$disturbance==scenarios$disturbance[z],"RecoveryRate"],results[results$disturbance==scenarios$disturbance[z],"collapsed"]))
  xSeq <- seq(min(xMat[,1]),max(xMat[,1]),by=0.1)
  polyMinMax <- cbind(c(xSeq,rev(xSeq)),c(sapply(xSeq,function(x){min(xMat[abs(xMat[which.min(abs(xMat[,1]-x)),1]-xMat[,1])<=0.025,2])}),sapply(rev(xSeq),function(x){max(xMat[abs(xMat[which.min(abs(xMat[,1]-x)),1]-xMat[,1])<=0.025,2])})))
  polygon(polyMinMax[,1],polyMinMax[,2],col=adjustcolor(dist_colours[z],transparency))}}))

plot(results$longMSY,results$longOcc,pch=21,bg=plotColours,xlab="Surplus production (25 years)",ylab="Patch occupancy (25 years)")
invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){polygon(findHull(df=results[results$disturbance==scenarios$disturbance[z],],x="longMSY",y="longOcc"),col=adjustcolor(dist_colours[z],transparency))}))

legend("bottomright",c("Uniform","Local, random","Local, extirpation"),pch=22,pt.bg=dist_colours,bty="n",title="Disturbance regime")
dev.off()


boxplot(longCV~stochastic+disturbance,data=results[results$spatial==scenarios$spatial[1] & results$temporal==scenarios$temporal[1],],col=dist_colours[results$disturbance],xlab="Stochasticity",ylab="Recovery rate (yr-1)")

tiff("Figures/Recovery along dispersal network and disturbance clines.tiff",compression="lzw",units="in",height=8,width=8,res=800)

matLayout <- matrix(0,nrow=16,ncol=16)
matLayout[1:8,1:8] <- 1
matLayout[1:8,9:16] <- 2
matLayout[9:16,1:8] <- 3
matLayout[9:16,9:16] <- 4
matLayout[5:7,6:8] <- 5
matLayout[5:7,14:16] <- 6
matLayout[13:15,6:8] <- 7
matLayout[13:15,14:16] <- 8

layout(matLayout)
par(mar=c(4,4,0.5,0.5))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i],]
  plot(subResults$dispersal,subResults$RecoveryRate,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery rate (yr-1)",xlab="Dispersal rate",ylim=1.0*range(c(0,results$RecoveryRate)),cex.lab=1.2)
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    #points(distResults$dispersal,distResults$recovery,bg=dist_colours[j],pch=21)
  }
  
  invisible(lapply(1:length(scenarios$disturbance),FUN=function(z){
    y_rng <- sapply(1:length(scenarios$dispersal),function(x){quantile(subResults$RecoveryRate[v_all.equal(subResults$dispersal,scenarios$dispersal[x]) & subResults$disturbance==scenarios$disturbance[z]],probs=c(0.25,0.75))});
    polygon(x=c(scenarios$dispersal,rev(scenarios$dispersal)),y=c(y_rng[1,],rev(y_rng[2,])),col=adjustcolor(dist_colours[z],transparency))}))
}
legend("right",c("Uniform","Local, random","Local, extirpation"),pch=22,pt.bg=dist_colours,bty="n",title="Disturbance regime",cex=1.2)
add2plot()
dev.off()

layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i],]
  vioplot(recovery~alpha,data=subResults,ylim=1.7*range(c(0,results$recovery)),drawRect=FALSE,col=c("orange","dodgerblue"))
  boxplot(recovery~alpha,data=subResults,ylim=1.7*range(c(0,results$recovery)),add=TRUE,col="grey30",boxwex=0.25,xaxt="n",yaxt="n",border="white",outpch=21,outcol="white",whisklty="solid",staplelwd=NA,outcex=1.5,outbg="grey30")
  abline(h=25,lwd=2,col="black",xpd=FALSE)
  mtext("Local productivity",1,line=2.5,cex=0.7)
  mtext("Recovery rate (generations)",2,line=2.5,cex=0.7)
}
add2plot()

# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
results$patchScen <- paste(results$alpha,"\u03B1","\n",results$beta,"\u03B2")
results$stochasticity <- as.character(results$stochastic)
p <- ggplot(data=results, aes(x=patchScen, y=RecoveryRate, fill=stochasticity))+ 
  geom_violin(trim=FALSE,scale = "width")+
  facet_wrap(~stochasticity)+
  geom_boxplot(width=0.1,colour="grey90",outlier.shape=NA,fill="grey50")+
  scale_fill_brewer(palette="Dark2") +
  labs(x="Patch variation",y="Recovery rate (yr-1)",fill="Stochasticity") +
  theme_minimal() +
  theme(legend.position="top",strip.text.x = element_blank())
p
ggsave("Figures/Recovery along variable local productivities.tiff",plot=p,compression="lzw",units="in",height=4,width=6,dpi=800)

layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery rate (generations)",xlab="Dispersal rate",ylim=1.5*range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery rate (generations)",xlab="Dispersal rate",ylim=1.5*range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[3] & 
                          results$spatial==scenarios$spatial[3] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$recovery,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Recovery rate (generations)",xlab="Dispersal rate",ylim=1.5*range(c(0,results$recovery)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$recovery,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY patterns
# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# MSY
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[2] & 
                          results$spatial==scenarios$spatial[2] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medMSY,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Sustainable production",xlab="Dispersal rate",ylim=1.5*range(c(0,results$medMSY)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medMSY,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# Capacity
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$med_capacity,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Carrying capacity",xlab="Dispersal rate",ylim=1.5*range(c(0,results$med_capacity)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$med_capacity,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# Compensation
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[1] & 
                          results$temporal==scenarios$temporal[1] & 
                          results$spatial==scenarios$spatial[1] & 
                          results$alpha==scenarios$alpha[1] & 
                          results$beta==scenarios$beta[1],]
  plot(subResults$dispersal,subResults$med_compensation,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Compensation",xlab="Dispersal rate",ylim=1.5*range(c(0,results$med_compensation)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$med_compensation,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()

# Occupancy
layout(matLayout)
par(mar=c(5,4,1,1))
for(i in 1:length(scenarios$network))
{
  subResults <- results[results$network==scenarios$network[i] & 
                          results$stochastic==scenarios$stochastic[2] & 
                          results$temporal==scenarios$temporal[2] & 
                          results$spatial==scenarios$spatial[2] & 
                          results$alpha==scenarios$alpha[2] & 
                          results$beta==scenarios$beta[2],]
  plot(subResults$dispersal,subResults$medOcc,col=0,lty=0,type="b",lwd=0,pch=0,ylab="Patch occupancy",xlab="Dispersal rate",ylim=1.70*range(c(0,results$medOcc)))
  for(j in 1:length(scenarios$disturbance)){
    distResults <- subResults[subResults$disturbance==scenarios$disturbance[j],]
    lines(distResults$dispersal,distResults$medOcc,col=dist_colours[j],type="b",pch=21,lwd=2,lty=1)
  }
}
add2plot()