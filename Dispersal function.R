
#recruitment <- rpois(25,5000)
dispersal <- function(maxDispersal,distDecay,dist_matrix,recruitment) # recruitment is named vector
{
  maxDispersal <- maxDispersal # dispersal rate
  distDecay <- distDecay # distance penalty for migration
  Npops <- ncol(dist_matrix)
  migrationNumer <- apply(dist_matrix,2,function(x){exp(-distDecay*x)}) # dispersal between sites, Anderson et al. 2015 EcoApps
  diag(migrationNumer) <- 0 # fish cannot migrate to their own population
  migration <- maxDispersal*(migrationNumer/colSums(migrationNumer))
  # rows correspond to emigrants away from patch
  # columns correspond to immigrants into patch
  dispersers <- migration*recruitment
  immigrants <- colSums(dispersers)
  emigrants <- rowSums(dispersers)
  return(list("immigrants"=immigrants,"emigrants"=emigrants,"dispersers"=dispersers))
}



