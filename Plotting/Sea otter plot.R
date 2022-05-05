library(ggplot2)
library(tidyr)

otters <- read.csv("Figures/SeaOtters_AKCA.csv",header=TRUE)
head(otters,2)
p1 <- ggplot(otters,aes(x=Year_Seq,y=Abundance,fill=Region,colour=Region)) +
  geom_ribbon(data=otters[otters$Method=='Estimated',],aes(x=Year_Seq,ymin=LCI,ymax=UCI,group=Region),colour=NA,fill='grey80',alpha=1) +
  geom_line(data=otters[otters$Method=='Estimated',],lwd=0.75) + 
  geom_hline(yintercept=3190,lty=2,lwd=0.75,col='black')+
  geom_point(data=otters[otters$Method=='Observed',],pch=21,col="black",size=0.75) +
  xlab("Years since initation of recovery") + ylab("Sea otter abundances") +
  theme_minimal() +
  annotate("text",label="California recovery target",x=10,y=15000,size=2.5)+
  geom_curve(x= 5, y = 14500,xend = 10, yend = 3500, color = 'black',arrow = arrow(length=unit(0.05,'inches'),type='closed'),lwd=0.5)+
  scale_fill_brewer(type="qual",palette=2) +
  scale_colour_brewer(type="qual",palette=2) +
  theme(legend.position="top")+
  geom_rect(aes(xmin=36,xmax=73,ymin=0,ymax=5000),col='black',lwd=0.25,fill=NA)
p1
sub_ott <- otters[otters$Region=='California' & otters$Year>=1983,]
ca <- ggplot(sub_ott,aes(x=Year_Seq,y=Abundance,fill=Region,colour=Region)) +
  geom_ribbon(data=sub_ott[sub_ott$Method=='Estimated',],aes(x=Year_Seq,ymin=LCI,ymax=UCI,group=Region),colour=NA,fill='grey80') +
  geom_line(data=sub_ott[sub_ott$Method=='Estimated',],lwd=0.2) + 
  geom_hline(yintercept=3190,lty=2,lwd=0.2,col='black')+
  geom_point(data=sub_ott[sub_ott$Method=='Observed',],pch=21,col="black",size=0.75) +
  theme_minimal() +
  xlab("") + ylab("") +
  annotate("text",label="(a)",x=38,y=2950,size=2)+
  scale_fill_manual(values="#d95f02") +
  scale_colour_manual(values="#d95f02") +
  theme(legend.position="none",panel.background = element_rect(fill=adjustcolor("white",0.90), colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text=element_text(size=5),axis.text=element_text())
ca
library(cowplot)
otter_inset <- ggdraw(p1) +
  draw_plot({
    ca},
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.65, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.60,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.35, 
    height = 0.25)
otter_inset
ggsave(plot=otter_inset,"Figures/Sea otter comparisons.jpeg",dpi=600,units='in',height=4.5,width=4.5)
