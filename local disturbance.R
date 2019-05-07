Disturbance <- function(metaPop,magnitude=0.5,DisType="uniform",N_p,prod,Npatches)
{
  # function to apply the disturbance regime to patchy populations
  metaPop = round(metaPop) # what is our total metapopulation
  totalLoss = round(metaPop*magnitude) # how much is our targeted loss across metapopulation
  N_p = round(N_p) # what is the current population size at the time of disturbance
  targetPatchMax <- min(which(cumsum(sort(N_p,decreasing=FALSE))>=totalLoss)) # calculate the minimum number of patches necessary to get target total losses
  targetPatchMin <- min(which(cumsum(sort(N_p,decreasing=TRUE))>=totalLoss)) # calculate the maximum number of patches necessary to get target total losses
  targetPatches <- round(mean(c(targetPatchMin,targetPatchMax))) # average number of matches necessary to get total loss
  animals_all <- data.frame("Animal"=1:sum(N_p),"Patch"=rep(1:Npatches,times=N_p)) # arrange animals as individuals susceptible to disturbance
  
  if(DisType=="uniform")
  {
    # randomly remove animals from populations, equal vulnerability across all patches
    target <- 0
    while(target<totalLoss)
    {
      total_p <- aggregate(Animal~Patch,data=animals_all,FUN=length)
      target <- nrow(animals_all)
    }
    harvest <- sample(1:nrow(animals_all),totalLoss)
    total_p <- animals_all[harvest,]
    mort_p <- aggregate(Animal~Patch,data=total_p,FUN=length)
    deaths_p <- rep(0,Npatches)
    deaths_p[mort_p$Patch] <- mort_p$Animal
  }
  
  if(DisType=="random")
  {
    # randomly remove animals from populations, equal vulnerability across selected patches
    target <- 0
    d_patches <- sample(1:Npatches,targetPatches)
    animals <- animals_all[animals_all$Patch%in%d_patches,]
    total_p <- aggregate(Animal~Patch,data=animals,FUN=length)
    target <- nrow(animals)
    while(target<totalLoss)
    {
      d_patches <- sample(1:Npatches,targetPatches)
      animals <- animals_all[animals_all$Patch%in%d_patches,]
      total_p <- aggregate(Animal~Patch,data=animals,FUN=length)
      target <- nrow(animals)
    }
    harvest <- sample(1:nrow(animals),totalLoss)
    total_p <- animals[harvest,]
    mort_p <- aggregate(Animal~Patch,data=total_p,FUN=length)
    deaths_p <- rep(0,Npatches)
    deaths_p[mort_p$Patch] <- mort_p$Animal
  }
  
  if(DisType=="random_patch")
  {
    # remove all animals from randomly selected patches to reach totalLoss
    target <- 0
    d_patches <- sample(1:Npatches,targetPatches)
    deaths_p <- rep(0,Npatches)
    deaths_p[d_patches] <- N_p[d_patches]
    target <- sum(deaths_p,na.rm=TRUE)
    while(target<totalLoss)
    {
      d_patches <- sample(1:Npatches,targetPatches)
      deaths_p <- rep(0,Npatches)
      deaths_p[d_patches] <- N_p[d_patches]
      target <- sum(deaths_p,na.rm=TRUE)
    }
  }
  
  if(DisType=="targeted")
  {
    # randomly remove animals from populations, vulnerability across selected patches depends on local productivity
    target <- 0
    d_patches <- sample(1:Npatches,targetPatches,prob=prod/sum(prod))
    animals <- animals_all[animals_all$Patch%in%d_patches,]
    total_p <- aggregate(Animal~Patch,data=animals,FUN=length)
    target <- nrow(animals)
    while(target<totalLoss)
    {
      d_patches <- sample(1:Npatches,targetPatches,prob=prod/sum(prod))
      animals <- animals_all[animals_all$Patch%in%d_patches,]
      total_p <- aggregate(Animal~Patch,data=animals,FUN=length)
      target <- nrow(animals)
    }
    harvest <- sample(1:nrow(animals),totalLoss)
    total_p <- animals[harvest,]
    mort_p <- aggregate(Animal~Patch,data=total_p,FUN=length)
    deaths_p <- rep(0,Npatches)
    deaths_p[mort_p$Patch] <- mort_p$Animal
  }
  # return the deaths by patch vector
  return(list("deaths_p"=deaths_p))
}