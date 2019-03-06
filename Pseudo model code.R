
source("Complex network.R")

distance_matrix <- distances(complexLandscape,v=V(complexLandscape),to=V(complexLandscape))

omega = 0.075
m = 0.15
alpha <- 20
beta <- 5e-3
cv <- 0.1
stock <- 0:1000
rec_mean <- alpha*stock*exp(-beta*stock)
rec.obs <- pmax(0,rnorm(length(stock),mean=rec_mean,sd=rec_mean*cv))

ricker <-function(Nadults){alpha*Nadults*exp(-beta*Nadults)}

plot(stock,rec_mean,type="l",ylim=c(0,max(rec.obs)))
points(stock,rec.obs,pch=21,bg="grey50")

# example dispersal outcome
dispersers_t <- dispersal(maxDispersal=omega,distDecay=m,dist_matrix,recruitment)
dispersers_t$immigrants

# example model
Nyears <- 100
Nstart <- 1000
Em <- rpois(Nyears,1)
Im <- rpois(Nyears,1)
Nadults <- rep(NA,Nyears)
Nrec <- rep(NA,Nyears)
# initialize at year 1
Nadults[1] <- Nstart
Nrec[1] <- ricker(Nadults[1])
pseudoSink <- sink <- source <- potentialRec <- rep(NA,Nyears)
for(II in 2:(Nyears+1))
{
  Nrec[II] <- ricker(Nadults[II-1])
  Nadults[II] <- Nrec[II] - Em[II] + Im[II]
  potentialRec[II] <- Nadults[II-1]-Im[II-1] # debate occurred.
  
  sourceRec[II] <- Nadults[II-1]-Im[II-1] # debate occurred.
  
  
  sink[II] <- (Nrec[II] < Nadults[II-1]) & (Em[II] < Im[II])
  source[II] <- (Nrec[II] > Nadults[II-1]) & (Em[II] > Im[II])
  pseudoSink[II] <- ((Nrec[II] < Nadults[II-1]) & (Nrec[II] > potentialRec[II])) & (Em[II] < Im[II])

}

plot(Nadults[1:Nyears],Nrec[2:(Nyears+1)])
abline(b=1,a=0,xpd=F,lty=2,col="red")

plot(sink,type="l")
sum(sink,na.rm=T)
sum(source,na.rm=T)
