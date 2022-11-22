results <- readRDS(file="Figures/clustering_outcomes.rds")
results$total <- 1
#as.numeric(table(paste(results$network_lab,results$disturb_lab,results$dispersal,results$patchScen,results$stochastic,results$spatial,results$temporal)))
surprises <- data.frame(aggregate(total~network_lab+disturb_lab+dispersal+patchScen+stochastic+spatial+temporal,data=results,FUN=sum))

results$cluster_numeric <- as.numeric(results$cluster_surprises)
results$patchScen <- gsub(" \n",",",results$patchScen)
# 
# HydeNet::inputCPT(cluster_surprises~network_lab+disturb_lab+dispersal_range+patchScen+stochastic_lab+space_var+temporal_var,
#                   data=results,
#                   factorLevels <- list(cluster_surprises=levels(results$cluster_surprises),
#                                        network_lab=levels(results$network_lab),
#                                        disturb_lab=levels(results$disturb_lab),
#                                        dispersal_range=levels(results$dispersal_range),
#                                        patchScen=levels(as.factor(results$patchScen)),
#                                        stochastic_lab=levels(results$stochastic_lab),
#                                        space_var=levels(results$space_var),
#                                        temporal_var=levels(results$temporal_var)),
#                   reduce=FALSE)
factorLevels <- list(cluster_surprises=levels(results$cluster_surprises),network_lab=levels(results$network_lab),disturb_lab=levels(results$disturb_lab),dispersal_range=levels(results$dispersal_range),patchScen=levels(as.factor(results$patchScen)),stochastic_lab=levels(results$stochastic_lab),space_var=levels(results$space_var),temporal_var=levels(results$temporal_var))


results_table <- addmargins(table(cluster_surprises=(results$cluster_surprises),
                                  network_lab=(results$network_lab),
                                  disturb_lab=(results$disturb_lab),
                                  dispersal_range=(results$dispersal_range),
                                  patchScen=(as.factor(results$patchScen)),
                                  stochastic_lab=(results$stochastic_lab),
                                  space_var=(results$space_var),
                                  temporal_var=(results$temporal_var)))
prob_res <- prob::probspace(results)
prob::Prob(prob_res,event=cluster_surprises=="Resilient",given=dispersal_range=="Low")
prob::Prob(prob_res,event=cluster_surprises=="Resilient",given=dispersal_range=="High")
prob::Prob(prob_res,event=cluster_surprises=="Resilient",given=dispersal_range=="Low" & network_lab=="Linear")
prob::Prob(prob_res,event=cluster_surprises=="Resilient",given=dispersal_range=="Low" & network_lab=="Grid")
prob::Prob(prob_res,event=cluster_surprises=="Slow recovery",given=dispersal_range=="Low" & network_lab=="Linear")
prob::Prob(prob_res,event=cluster_surprises=="Slow recovery",given=dispersal_range=="Low" & network_lab=="Grid")

olr <- MASS::polr(cluster_surprises~network_lab+disturb_lab+dispersal_range+patchScen+stochastic_lab+space_var+temporal_var,data=results,Hess=TRUE)
summary(olr)

library("effects")
Effect(focal.predictors = "network_lab",olr)
Effect(focal.predictors=c("network_lab","disturb_lab","dispersal_range","patchScen","stochastic_lab","space_var","temporal_var"),olr)
plot(Effect(focal.predictors = "network_lab",olr))
plot(Effect(focal.predictors = c("network_lab", "disturb_lab",'dispersal_range'),olr))
plot(Effect(focal.predictors = c("patchScen","network_lab", "stochastic_lab"),olr))
plot(Effect(focal.predictors = c("patchScen","network_lab", "disturb_lab"),olr),lattice=list(key.args=list(cex=0,cex.title=0.5)))
plot(Effect(focal.predictors = c("patchScen","network_lab", "disturb_lab"),olr,factor.names=FALSE),lattice=list(key.args=list(cex=2,cex.title=1)),axes=list(grid=TRUE,x=list(rotate=30)))

gg_resp <- ggeffects::ggpredict(olr, terms = c("network_lab","disturb_lab","dispersal_range","patchScen"),
                                condition =c("stochastic_lab"="High variance",
                                             "space_var"="High spatial",
                                             "temporal_var"="High temporal"))

gg_sub <- gg_resp[gg_resp$response.level=="Resilient",]
attr(gg_sub, "title") <- "Predicted probabilities of resilient recoveries"
attr(gg_sub, "x.title") <- "Habitat network"
attr(gg_sub, "y.title") <- "Predicted probabilities"
attr(gg_sub, "legend.title")= "Disturbance regime"
plot(gg_sub,facet=TRUE,ci=FALSE,colors="system", connect.lines = TRUE)
library(ggplot2)
ggplot(gg_resp[gg_resp$response.level=="Resilient",],aes(x=x,y=predicted,colour=group))+
  geom_point() +
  geom_line() +
  facet_grid(row=vars(panel),cols=vars(facet)) +
  theme_minimal()+
  scale_color_brewer()
