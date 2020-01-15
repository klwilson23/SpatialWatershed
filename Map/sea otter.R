library(sf)
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
---------------------------##SEA OTTER FIGURE###-----------------
#########################################################################

# Set working directory
getwd()
setwd("D:/Lauren files/Sea Otter/Shapefiles")


# read in california shapefile and plot california map
california<-readOGR("california.shp")
plot(california)

california_map <- ggplot() + 
  geom_polygon(data = california, aes(x=long, y = lat, group = group), colour = "grey70", fill = "grey96") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(california_map)



# read in polygons shapefile and plot polygon map
polygons<-readOGR("CA_polygons.shp")
plot(polygons)

polygons

polygons_map <- ggplot() + 
  geom_polygon(data = polygons, aes(x=long, y = lat, group = group, colour = Pop) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(polygons_map)

?group

# combine california and polygons shapefiles into one plote
combine <- ggplot() + 
  geom_polygon(data = california, aes(x=long, y = lat, group = group), colour = "grey70", fill = "grey96") +
  geom_polygon(data = polygons, aes(x=long, y = lat, group = group), colour = "pink", fill = "pink") +
  labs(title = "Sea Otter Movement", x = "Longitude", y="Latitude") +
  coord_fixed(ratio = 1, xlim = c(-300000,90000), ylim = c(00000,-460000), expand = TRUE, clip = "on")
plot(combine)


# read sea otter density data into R
density <- read.csv(file = "poptrends.csv", header = TRUE)
colnames(density)
density


density_map <- ggplot() + 
  geom_polygon(data = density, aes(x=long, y = lat, group = group), colour = Density) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(density_map)


lat <- pts$lat
long <- pts$long
population <- pts$population

map.points <- combine + geom_point(data = pts, aes(x = long, y = lat, colour = population), pch=16, size=3.5, alpha=I(0.9)) +
  scale_color_manual(breaks = c("northern", "southern"), values=c("royalblue3", "orange2"))
plot(map.points)

