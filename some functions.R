Corner_text <- function(text, location="topright",...) #function to write text to the corner of plots
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}

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

dvSpaceTime <- function(mnSig,lastDV,rhoTime,rhoSpace,distMatrix)
{
  # function to calculate recruitment deviations correlated through space and time
  ar1 <- rhoTime*lastDV
  dvST <- ar1 + rmvnorm(1,mean=rep(0,ncol(distMatrix)),sigma=(mnSig*(1-rhoTime^2)*(exp(-rhoSpace*distMatrix))))
  return(dvST)
}

SRfn <- function(theta){
  a.hat <- theta[1]
  b.hat <- exp(theta[2])
  sd.hat <- exp(theta[3])
  rec.mean <- (a.hat*spawnRec$spawners)/(1+((a.hat-1)/b.hat)*spawnRec$spawners)
  nll <- -1*sum(dlnorm(spawnRec$recruits,meanlog=log(rec.mean),sdlog=sd.hat,log=TRUE)*spawnRec$weights)
  penalty1 <- -dnorm(a.hat,alpha,3*alpha,log=TRUE)
  penalty2 <- -dnorm(b.hat,metaK,3*metaK,log=TRUE)
  jnll <- nll +penalty1+penalty2
  return(nll)
}
