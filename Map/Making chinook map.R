library(sf)
library(tidyverse)
library(rgdal)
library(ggplot2)
library(tmap)
library(tmaptools)
library(leaflet)
library(dplyr)

library("maps")
library("maptools")
library("raster")
library("sp")
library("foreign")
library("grid")
library("mapproj")
library("rgeos")
library("ggmap")
library("lattice")
library("png")

library(car)
library(gridExtra)
library(Kendall)
library(lmtest)
library(randtests)
library(zyp)

#########################################################################
----------------##SNAKE RIVER FIGURE###--------------------
#########################################################################

# Set working directory

# read in snake river states shapefile and plot states map
ch_srf<-readOGR("C:/Users/kylel/Google Drive/Project_MooreLabPaper/Literature/Case studies/B. Salmon/Figures/Chinook/CKSRF.shp")

ch_srs<-readOGR("C:/Users/kylel/Google Drive/Project_MooreLabPaper/Literature/Case studies/B. Salmon/Figures/Chinook/CKSRS.shp")

plot(ch_srf)
str(ch_srf)
state_map <- ggplot() + 
  geom_polygon(data = ch_srf, aes(x=long, y = lat, group = group), colour = "grey70", fill = "grey96") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(state_map)


combine <- ggplot() + 
  geom_polygon(data = ch_srf, aes(x=long, y = lat, group = STATUS), colour = "grey70", fill = "lightblue") +
  geom_polygon(data = ch_srs, aes(x=long, y = lat, group = STATUS), colour = "lightblue", fill = "lightblue") +
  labs(title = "Snake River Chinook salmon", x = "Longitude", y="Latitude") +
  coord_fixed(ratio = 1, expand = TRUE, clip = "on")
plot(combine)
