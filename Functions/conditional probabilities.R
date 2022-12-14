library(ggplot2)
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
#results$network_lab <- relevel(results$network_lab,"Grid")
scenarios <- list(cluster_surprises=levels(results$cluster_surprises),network_lab=levels(results$network_lab),disturb_lab=levels(results$disturb_lab),dispersal_range=levels(results$dispersal_range),patchScen=levels(as.factor(results$patchScen)),stochastic_lab=levels(results$stochastic_lab),space_var=levels(results$space_var),temporal_var=levels(results$temporal_var))
scenarios <- expand.grid(scenarios)
cond_prob <- tidyr::pivot_longer(scenarios,cols=2:8,names_to='scenario',values_to="condition")
cond_prob$prob <- NA
cond_prob$inv_prob <- NA
cond_prob <- cond_prob[!duplicated(cond_prob),]

for(i in 1:nrow(cond_prob))
{
  cond <- as.character(cond_prob$condition[i])
  label <- intersect(colnames(prob_res),cond_prob$scenario[i])
  label_names <- c("probs","cluster_surprises",label)
  sub_prob <- prob_res[,label_names]
  colnames(sub_prob) <- c("probs","cluster_surprises","scenario")
  cond_prob$prob[i] <-sum(sub_prob$probs[sub_prob$cluster_surprises==cond_prob$cluster_surprises[i] & sub_prob$scenario==cond])/(sum(sub_prob$scenario==cond)/nrow(sub_prob))
  cond_prob$inv_prob[i] <-sum(sub_prob$probs[sub_prob$cluster_surprises!=cond_prob$cluster_surprises[i] & sub_prob$scenario==cond])/(sum(sub_prob$scenario==cond)/nrow(sub_prob))
}
#cond_prob$condition <- as.character(cond_prob$condition)
ggplot(cond_prob,aes(x=condition,y=prob))+
  geom_point(data=cond_prob,aes(x=condition,y=prob),pch=21,fill="purple",colour="black") +
  geom_line(data=cond_prob,aes(x=condition,y=prob,group=interaction(cluster_surprises,scenario)),colour="purple") +
  facet_grid(row=vars(cluster_surprises),cols=vars(scenario),scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35,hjust = 1))

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

scenarios <- list(network_lab=levels(results$network_lab),disturb_lab=levels(results$disturb_lab),dispersal_range=levels(results$dispersal_range),patchScen=levels(as.factor(results$patchScen)),stochastic_lab=levels(results$stochastic_lab),space_var=levels(results$space_var),temporal_var=levels(results$temporal_var))
scenarios <- expand.grid(scenarios)
scenarios <- cbind(scenarios,predict(olr,newdata=scenarios,type="probs"))
scenarios <- tidyr::pivot_longer(scenarios,cols=levels(results$cluster_surprises),names_to="cluster_surprises",values_to="probs")
scenarios$variance <- paste(scenarios$stochastic_lab,scenarios$space_var,scenarios$temporal_var,sep=", ")
sub_scen <- scenarios[scenarios$variance=="High variance, High spatial, High temporal" | 
                        scenarios$variance=="Low variance, Low spatial, Low temporal",]
resilient <- sub_scen[sub_scen$cluster_surprises=="Resilient",]


make_plot <- function(data,xxx)
{
  plot_data <- data[sub_scen$cluster_surprises==xxx,]
  if(xxx=="Resilient")
  {
    plot_label <- paste("Probability of",tolower(xxx),"metapopulation recoveries",sep=" ")
  }else{
    plot_label <- paste("Probability of",tolower(xxx),"in metapopulation recoveries",sep=" ")
  }
  ggplot(data=plot_data,aes(x=network_lab,y=probs,fill=disturb_lab,shape=variance))+
    geom_line(data=plot_data,aes(x=network_lab,y=probs,group=interaction(variance,disturb_lab,patchScen,dispersal_range))) +
    geom_point(data=plot_data,aes(x=network_lab,y=probs,fill=disturb_lab,shape=variance)) +
    facet_grid(row=vars(patchScen),cols=vars(dispersal_range),scales="fixed") +
    theme_minimal()+
    #scale_color_brewer(name="Disturbance regime")+
    scale_fill_brewer(name="Disturbance regime",type="qual",palette=2)+
    scale_shape_manual(name="Variance scenario",values=c(21,22),labels=c("High","Low")) +
    labs(title=plot_label)+
    xlab("Habitat network")+ylab("Probability of outcome")+
    theme(legend.position="top",panel.spacing = unit(2, "lines"))+
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  ggsave(paste("Figures/",plot_label,".jpeg",sep=""),dpi=600,units='in',height=8,width=8)
}
make_plot(sub_scen,"Resilient")
make_plot(sub_scen,"Slow recovery")
make_plot(sub_scen,"Hidden collapses")
make_plot(sub_scen,"Lost capacity")
make_plot(sub_scen,"Critical risk")
