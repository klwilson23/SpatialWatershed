patch_variance <- function(alpha_heterogeneity=F,cap_heterogeneity=F,Npatches=25,alpha_p,k_p)
{
  # replace internal function with Colin's RMD file dynamics
  if(alpha_heterogeneity)
  {
    alpha_M <- sum(alpha_p)
    sample_prop <- runif(Npatches,0,1)
    alpha_p <- alpha_M*(exp(sample_prop)/sum(exp(sample_prop)))
    beta_p <- -log(alpha_p)/k_p
  }
  
  if(cap_heterogeneity)
  {
    sample_prop <- runif(Npatches,0,1)
    k_p <- metaK*(exp(sample_prop)/sum(exp(sample_prop)))
    beta_p <- -log(alpha_p)/k_p
  }
  return(list("alpha_p"=alpha_p,"beta_p"=beta_p,"k_p"=k_p))
}

# ignore this below

function()
{
  metaRec <- rowSums(rec_mean)
  metaStock <- stock*Npatches
  
  # test to see if we recover the starting parameters
  metaRicker <- lm(log(metaRec/metaStock)~metaStock)
  
  SR_df <- data.frame("stock"=metaStock,"recruits"=metaRec)
  
  fit_nls <- nls(recruits~alpha*stock*exp(-(log(alpha)/K)*stock),data=SR_df,start=list("alpha"=alpha,"K"=metaK))
  
  meta_alpha <- exp(coef(metaRicker)[1])
  meta_beta <- coef(metaRicker)[2]
  -log(meta_alpha)/meta_beta
  
  plot(metaStock,metaRec,type="l",lwd=2,col="dodgerblue")
  invisible(sapply(1:Npatches,function(x){
    lines(stock,rec_mean[,x],type="l",col="orange")
  }))
  
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
  
}