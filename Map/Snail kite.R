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
----------------##SNAIL KITE FIGURE###--------------------
#########################################################################

# Set working directory
getwd()
setwd("D:/Lauren files/Snail Kite/Shapefiles")


# read in florida shapefile and plot florida map
florida<-readOGR("florida_proj.shp")
plot(florida)

florida_map <- ggplot() + 
  geom_polygon(data = florida, aes(x=long, y = lat, group = group), colour = "grey39", fill = "grey39") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(florida_map)


# read in lake okeecheeba shapefile and plot lakeokee map
lakeokee<-readOGR("okeeproj.shp")
plot(lakeokee)

lakeokee_map <- ggplot() + 
  geom_polygon(data = lakeokee, aes(x=long, y = lat, group = group), colour = "grey39", fill = "white") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(lakeokee_map)


# read in lake jackson shapefile and plot lakejack map
lakejack<-readOGR("jacksonproj.shp")
plot(lakejack)

lakejack_map <- ggplot() + 
  geom_polygon(data = lakejack, aes(x=long, y = lat, group = group), colour = "grey39", fill = "white") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(lakejack_map)



# combine florida, lakeokee, lakejack shapefiles into one plot
combine <- ggplot() + 
  geom_polygon(data = florida, aes(x=long, y = lat, group = group), colour = "grey70", fill = "grey96") +
  geom_polygon(data = lakeokee, aes(x=long, y = lat, group = group), colour = "slategray2", fill = "slategray2") +
  geom_polygon(data = lakejack, aes(x=long, y = lat, group = group), colour = "slategray2", fill = "slategray2") +
  labs(title = "Snail Kite Movement", x = "Longitude", y="Latitude") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(combine)



# read point data into R
pts <- read.csv(file = "points.csv", header = TRUE)
colnames(pts)
pts

lat <- pts$lat
long <- pts$long
population <- pts$population

map.points <- combine + geom_point(data = pts, aes(x = long, y = lat, colour = population), pch=16, size=3.5, alpha=I(0.9)) +
  scale_color_manual(breaks = c("northern", "southern"), values=c("royalblue3", "orange2"))
plot(map.points)


#######

utmcoor<-SpatialPoints(cbind(combine$Easting, combine$Northing), proj4string=CRS("+proj=utm +zone=17"))
longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))


