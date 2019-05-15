source("make networks.R")
source("Dispersal function.R")
source("patch_variance.R")
source("local disturbance.R")
source("some functions.R")
source("popDynFn.R")

######################
## Differentiation ###
######################

yieldRoot<-function(alpha,beta,Fcur,model)
{
  Ucur <- 1-exp(-Fcur)
  if(model=="Beverton-Holt")
  {
    # population model is R ~ (CR * S)/(1+(CR-1)/K * S): reparameterized Beverton-Holt
    adults <- beta*exp(-Fcur)
    recruits <- (alpha*adults)/(1+((alpha-1)/beta)*adults)
    yield <- recruits-adults
  }
  if(model=="Ricker")
  {
    # population model is R ~ CR * S * e(-log(CR)/K * S): reparameterized Ricker
    adults <- (-log(alpha)/beta)*exp(-Fcur)
    recruits <- alpha*adults*exp(beta*adults)
    yield <- recruits-adults
  }
  return(list(recruits=recruits,adults=adults,yield=yield))
}

fun_yield <- function(Fcur,alpha,beta,delta,model)
{
  y1 <- yieldRoot(alpha=alpha,beta=beta,Fcur=Fcur-delta/2,model=model)$yield
  y2 <- yieldRoot(alpha=alpha,beta=beta,Fcur=Fcur+delta/2,model=model)$yield
  approx.gradient <- (y2-y1)/delta
  return(approx.gradient)
}

###############################################
## uniroot the differentiation for Fmsy #######
###############################################
findMSY <- function(alpha,beta,Npatches,model){
  F_msy <- NULL
  for(i in 1:Npatches){
    a_p <- alpha[i]
    b_p <- beta[i]
    F_msy[i] <- uniroot(fun_yield, interval=c(0,1),extendInt="yes",alpha=a_p,beta=b_p, delta=0.0001,model=model)$root 
  }
  MSY <- sapply(1:Npatches,function(x){yieldRoot(alpha=alpha[x],beta=beta[x],Fcur=F_msy[[x]][1],model=model)})
  return(list("F_msy"=F_msy,"MSY"=MSY))
}

Npatches <- 5
alpha <- runif(Npatches,1.2,3)
beta <- runif(Npatches,50,150)
model <- "Beverton-Holt"
findMSY(alpha=alpha,beta=beta,Npatches=Npatches,model=model)


curve((alpha[1]*x)/(1+((alpha[1]-1)/beta[1])*x)-x,from=0,to=max(beta))
abline(v=findMSY(alpha=alpha,beta=beta,Npatches=Npatches,model=model)$MSY["adults",1])
