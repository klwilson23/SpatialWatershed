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
getwd()
setwd("D:/Lauren files/Snake River/Shapefiles")


# read in snake river states shapefile and plot states map
states<-readOGR("states_proj.shp")
plot(states)

state_map <- ggplot() + 
  geom_polygon(data = states, aes(x=long, y = lat, group = group), colour = "grey70", fill = "grey96") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(state_map)


# read in snake river shapefile and plot snake river map
snake <-readOGR("snake_proj.shp")
plot(snake)

snake_map <- ggplot() + 
  geom_polygon(data = snake, aes(x=long, y = lat, group = group), colour = "lightblue", fill = "lightblue") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(snake_map)


# read in salmon river shapefile and plot salmon river map
salmon <-readOGR("salmon_proj.shp")
plot(salmon)

salmon_map <- ggplot() + 
  geom_polygon(data = salmon, aes(x=long, y = lat, group = group), colour = "lightblue", fill = "lightblue") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(salmon_map)


# read in columbia river shapefile and plot columbia river map
columbia <-readOGR("columbia_proj.shp")
plot(columbia)

columbia_map <- ggplot() + 
  geom_polygon(data = columbia, aes(x=long, y = lat, group = group), colour = "lightblue", fill = "lightblue") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
plot(columbia_map)

columbia


# read dam data into R
dams <- read.csv(file = "dams.csv", header = TRUE)
colnames(dams)
dams

lat <- dams$Latitude
long <- dams$Longitude
dam <- dams$dam

dams_map <- ggplot() + 
  geom_point(data = dams, aes(x=long, y=lat),size=5)
plot(dams_map)

####### convert dam data from lat and lon to UTM #########

head(dams, n=5)
sapply(dams, class)

# make the UTM cols spatial (X/Easting/Lon, Y/Northig/Lat
dams.SP <- st_as_sf(dams, coords = c("Longitude", "Latitude"), crs = 4326)

# transform to UTM
dams.SP <- st_transform(x = dams.SP, crs = 32612)

#get coordinates and add back to dataframe
dams.SP$utm_E <- st_coordinates(dams.SP)[,1]
dams.SP$utm_N <- st_coordinates(dams.SP)[,2]

# now switch back to Lat - Long
dams.SP <- st_transform(x = dams.SP, crs = 4326)

# add coordinates to dataframe
dams.SP$lon <- st_coordinates(dams.SP)[,1]
dams.SP$lat <- st_coordinates(dams.SP)[,2]

# coerce back to data frame
dams.SP<-st_set_geometry(dams.SP, NULL)

head(dams.SP, n=5)

# map dam points
dams_map <- ggplot() + 
  geom_point(data = dams.SP, aes(x = utm_E, y = utm_N), size=5)
plot(dams_map)

######END#######

####### combine states and rivers into one plot ####
combine <- ggplot() + 
  geom_polygon(data = states, aes(x=long, y = lat, group = group), colour = "grey70", fill = "grey96") +
  geom_polygon(data = salmon, aes(x=long, y = lat, group = group), colour = "lightblue", fill = "lightblue") +
  geom_polygon(data = snake, aes(x=long, y = lat, group = group), colour = "lightblue", fill = "lightblue") +
  geom_polygon(data = columbia, aes(x=long, y = lat, group = group), colour = "lightblue", fill = "lightblue") +
  geom_point(data = dams.SP, aes(x = utm_E, y = utm_N),  pch=15, size=2) +
  labs(title = "Snake River Salmon Movement", x = "Longitude", y="Latitude") +
  coord_fixed(ratio = 1, xlim = c(-600000,550000), ylim = c(4600000, 5500000), expand = TRUE, clip = "on")
plot(combine)


