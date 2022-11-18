library(mvtnorm)
library(vioplot)
library(ggplot2)
library(ggalluvial)
library("alphahull")
library(wesanderson)
library(ggpubr)
library(cluster)
library(fpc)
library(clValid)

wesAnderson <- "Moonrise1"

bootstrap_results <- readRDS("Simulations/bootstrap_results2022-11-17.rds")
bootstrap_results$StateShift <- ((1-bootstrap_results$recovered)+(1-bootstrap_results$longOcc))/2
bootstrap_results$RecoveryRate <- 1-bootstrap_results$recovery/NyrsPost
bootstrap_results$collapsed <- 1-bootstrap_results$recovered

bootstrap_results <- tidyr::pivot_longer(bootstrap_results,cols=colnames(bootstrap_results)[-c(1:8,35)],names_to="Metric",values_to="Value")
sub_results <- bootstrap_results[bootstrap_results$Metric%in%c("RecoveryRate","collapsed","metaAbund","short_surp","med_surp","long_surp","shortOcc","medOcc","longOcc"),]

sub_results$Metric <- factor(sub_results$Metric,levels=c("RecoveryRate","collapsed","metaAbund","short_surp","med_surp","long_surp","shortOcc","medOcc","longOcc"),labels=c("Recovery rate","Prop. of sims w/ no recovery","Metapopulation Abundance","Short-term production","Medium-term production","Long-term production","Short-term occupancy","Medium-term occupancy","Long-term occupancy"))

ggplot(sub_results,aes(x=Nboots,y=Value,col=Nboots)) +
  geom_point()+
  geom_line()+
  facet_wrap(~Metric) +
  xlab("Number of bootstrap iterations") + ylab("Value for metric")
