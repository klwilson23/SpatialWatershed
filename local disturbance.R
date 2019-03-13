Disturbance <- function(metaPop,magnitude=0.5,DisType="uniform",N_p,prod=alpha_p)
{
  metaPop = round(metaPop)
  totalLoss = round(metaPop*magnitude)
  N_p = round(N_p)
  targetPatchMax <- min(which(cumsum(sort(N_p,decreasing=FALSE))>=totalLoss)) # calculate the minimum number of patches necessary to get target total losses
  targetPatchMin <- min(which(cumsum(sort(N_p,decreasing=TRUE))>=totalLoss)) # calculate the maximum number of patches necessary to get target total losses
  targetPatch <- round(mean(c(targetPatchMin,targetPatchMax)))
  animals_all <- data.frame("Animal"=1:sum(N_p),"Patch"=rep(1:Npatches,times=N_p))
  
  if(DisType=="uniform")
  {
    # number of animals removed, random by patches
    target <- 0
    while(target<totalLoss)
    {
      total_p <- aggregate(Animal~Patch,data=animals_all,FUN=length)
      target <- nrow(animals_all)
    }
    harvest <- sample(1:nrow(animals_all),totalLoss)
    survivors <- animals_all[-harvest,]
    surv_p <- aggregate(Animal~Patch,data=survivors,FUN=length)
    deaths_p <- rep(0,Npatches)
    deaths_p[surv_p$Patch] <- (N_p[surv_p$Patch]-surv_p$Animal)
  }
  
  if(DisType=="random")
  {
    # same number of animals removed, patches are randomly chosen
    target <- 0
    d_patches <- sample(1:Npatches,targetPatches)
    animals <- animals_all[animals_all$Patch%in%d_patches,]
    while(target<totalLoss)
    {
      d_patches <- sample(1:Npatches,targetPatches)
      animals <- animals_all[animals_all$Patch%in%d_patches,]
      total_p <- aggregate(Animal~Patch,data=animals,FUN=length)
      target <- nrow(animals)
    }
    harvest <- sample(1:nrow(animals),totalLoss)
    survivors <- animals[-harvest,]
    surv_p <- aggregate(Animal~Patch,data=survivors,FUN=length)
    deaths_p <- rep(0,Npatches)
    deaths_p[surv_p$Patch] <- (N_p[surv_p$Patch]-surv_p$Animal)
  }
  
  if(DisType=="random_patch")
  {
    # randomized: half of the patches go extinct
    d_patches <- sample(1:Npatches,targetPatches)
    deaths_p <- rep(0,Npatches)
    deaths_p[d_patches] <- N_p[d_patches]
  }
  
  if(DisType=="targeted")
  {
    # same number of animals removed, with more productive patches most likely targeted
    target <- 0
    d_patches <- sample(1:Npatches,targetPatches,prob=exp(prod)/sum(exp(prod)))
    animals <- animals_all[animals_all$Patch%in%d_patches,]
    while(target<totalLoss)
    {
      d_patches <- sample(1:Npatches,targetPatches,prob=exp(prod)/sum(exp(prod)))
      animals <- animals_all[animals_all$Patch%in%d_patches,]
      total_p <- aggregate(Animal~Patch,data=animals,FUN=length)
      target <- nrow(animals)
    }
    harvest <- sample(1:nrow(animals),totalLoss)
    survivors <- animals[-harvest,]
    surv_p <- aggregate(Animal~Patch,data=survivors,FUN=length)
    deaths_p <- rep(0,Npatches)
    deaths_p[surv_p$Patch] <- (N_p[surv_p$Patch]-surv_p$Animal)
  }
  return(list("deaths_p"=deaths_p))
}