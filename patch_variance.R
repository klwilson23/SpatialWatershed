patch_variance <- function(model="Beverton-Holt",alpha_heterogeneity=F,cap_heterogeneity=F,Npatches=16,alpha,alpha_p,metaK,k_p,magnitude_of_decline)
{
  if(alpha_heterogeneity & model=="Ricker")
  {
    alpha_M <- alpha
    alpha_p <- rnorm(Npatches,0,1)
    alpha_p <- (alpha_p-mean(alpha_p))/sd(alpha_p)
    alpha_p <- alpha_M + alpha_p 
    while(any(alpha_p<1.0)){
      alpha_p <- rnorm(Npatches,0,1)
      alpha_p <- (alpha_p-mean(alpha_p))/sd(alpha_p)
      alpha_p <- alpha_M + alpha_p 
    }
    beta_p <- -log(alpha_p)/k_p
  }
  
  if(cap_heterogeneity & model=="Ricker")
  {
    sample_prop <- runif(Npatches,0,1)
    k_p <- as.vector(rmultinom(1,metaK,prob=(exp(sample_prop)/sum(exp(sample_prop)))))
    while(any(k_p >= metaK*(1-magnitude_of_decline))){
      k_p <- as.vector(rmultinom(1,metaK,prob=(exp(sample_prop)/sum(exp(sample_prop)))))
    }
    beta_p <- -log(alpha_p)/k_p
  }
  
  if(alpha_heterogeneity & model=="Beverton-Holt")
  {
    alpha_M <- alpha
    alpha_p <- rnorm(Npatches,0,1)
    alpha_p <- (alpha_p-mean(alpha_p))/sd(alpha_p)
    alpha_p <- alpha_M + alpha_p 
    while(any(alpha_p<1.0)){
      alpha_p <- rnorm(Npatches,0,1)
      alpha_p <- (alpha_p-mean(alpha_p))/sd(alpha_p)
      alpha_p <- alpha_M + alpha_p 
    }
    beta_p <- k_p
  }
  
  if(cap_heterogeneity & model=="Beverton-Holt")
  {
    sample_prop <- runif(Npatches,0,1)
    k_p <- as.vector(rmultinom(1,metaK,prob=(exp(sample_prop)/sum(exp(sample_prop)))))
    
    while(any(k_p >= metaK*(1-magnitude_of_decline))){
      k_p <- as.vector(rmultinom(1,metaK,prob=(exp(sample_prop)/sum(exp(sample_prop)))))
    }
    
    beta_p <- k_p
  }
  
  return(list("alpha_p"=alpha_p,"beta_p"=beta_p,"k_p"=k_p))
}

# ignore this below

function()
{
  stock <- 1:(2*max(k_p))
  rec_mean <- sapply(1:Npatches,function(x){(alpha_p[x]*stock)/(1+((alpha_p[x]-1)/beta_p[x])*stock)})

  metaRec <- rowSums(rec_mean)
  metaStock <- stock*Npatches
  
  SR_df <- data.frame("spawners"=metaStock,"recruits"=metaRec)
  plot(metaStock,metaRec,type="l",lwd=2,col="dodgerblue")
  invisible(sapply(1:Npatches,function(x){
    lines(stock,rec_mean[,x],type="l",col="orange")
  }))
  
  # test to see if we recover the starting parameters
  metaRicker <- lm(log(metaRec/metaStock)~metaStock)
  
  fit_nls <- nls(recruits~alpha*stock*exp(-(log(alpha)/K)*stock),data=SR_df,start=list("alpha"=alpha,"K"=metaK))
  
  meta_alpha <- exp(coef(metaRicker)[1])
  meta_beta <- coef(metaRicker)[2]
  -log(meta_alpha)/meta_beta
  
  matplot(stock,rec_mean)
  
  SR_df <- data.frame("stock"=metaStock,"recruits"=metaRec)
  fit_nls <- nls(recruits~alpha*stock*exp(-(log(alpha)/K)*stock),data=SR_df,start=list("alpha"=alpha,"K"=metaK),algorithm="port")
  
  ricker_SR <- function(l_alpha,l_K,l_sd){
    sd <- exp(l_sd)
    K <- exp(l_K)
    alpha <- exp(l_alpha)
    rec_mean <- alpha*SR_df$stock*exp(-(log(alpha)/K)*SR_df$stock)
    nll <- -(sum(dnorm(SR_df$recruits,mean=rec_mean,sd=sd,log=T),na.rm=T))
    return("nll"=nll)
  }
  theta <- list("l_alpha"=log(20),"l_K"=log(25000),"l_sd"=log(1000))
  mle_fit <- mle(ricker_SR,start=theta,method="SANN",nobs=nrow(SR_df))
  exp(mle_fit@coef)
  meta_alpha <- exp(coef(metaRicker)[1])
  meta_beta <- coef(metaRicker)[2]
  
  print(c(meta_alpha,meta_beta,-log(meta_alpha)/meta_beta))
  summary(fit_nls)
  
  #SR_df <- data.frame("stock"=SR_df$stock/max(SR_df$stock),"recruits"=SR_df$recruits/max(SR_df$recruits))
  SRfit <- nls(recruits~(a.hat*spawners)/(1+((a.hat-1)/b.hat)*spawners),data=SR_df,start=list("a.hat"=3,"b.hat"=1000),lower=c(0,0),upper=c(Inf,Inf),algorithm="port",control=list(tol=1e-12))
  alphaHat <- coef(SRfit)["a.hat"]
  metaK_hat <- coef(SRfit)["b.hat"]
  
  plot(recruits~spawners,data=SR_df,type="l",lwd=2,col="dodgerblue")
  curve((alphaHat*x)/(1+((alphaHat-1)/metaK_hat)*x),add=T,lwd=2,col="orange")
  curve((alpha*x)/(1+((alpha-1)/metaK)*x),add=T,lwd=2,col="red")
  
}