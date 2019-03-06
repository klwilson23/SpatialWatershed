
source("Complex network.R")
source("Dispersal function.R")

distance_matrix <- distances(complexLandscape,v=V(complexLandscape),to=V(complexLandscape))

ricker <-function(Nadults){alpha*Nadults*exp(beta*Nadults)}

Npatches <- ncol(distance_matrix)

# leading parameters
omega <- 0.075 # proportion of animals in patch that move
m <- 0.15 # distance decay function: could set at some proportion of max distance
# adult stock-juvenile recruitment traits
alpha <- 20
metaK <- 25000
beta <- -log(alpha)/metaK
cv <- 0.1
alpha_heterogeneity <- FALSE
beta_heterogeneity <- TRUE

alpha_p <- rep(alpha,Npatches)
k_p <- (metaK/Npatches)
beta_p <- rep(-log(alpha)/k_p,Npatches)

if(alpha_heterogeneity)
{
  alpha_M <- sum(alpha_p)
  sample_prop <- runif(Npatches,0,1)
  alpha_p <- alpha_M*(exp(sample_prop)/sum(exp(sample_prop)))
}

if(beta_heterogeneity)
{
  sample_prop <- runif(Npatches,0,1)
  k_p <- metaK*(exp(sample_prop)/sum(exp(sample_prop)))
  beta_p <- -log(alpha_p)/k_p
}

k_p <- -log(alpha_p)/beta_p
sum(k_p)

stock <- 0:(round(max(k_p))+1)
rec_mean <- sapply(1:Npatches,function(x){alpha_p[x]*stock*exp(beta_p[x]*stock)})

dim(rec_mean)

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

rec.obs <- pmax(0,rnorm(length(stock),mean=rec_mean,sd=rec_mean*cv))
