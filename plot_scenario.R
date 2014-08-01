# Plot scenario.R
################################
# Date: 31/07/2014
# Author: CNM 
# Project: TFM
# Description: Script to plot scenarios and bym, crossval, max, mean, wmean 
#              results.


# Clean workspace & set wd
rm(list=ls())
# setwd("~/GitHub/tfm")
filesDir <- "~/Dropbox/Claudio-Migue"
workDir <- "~/Documents/docs/projects/tfm"
setwd(workDir)

# Libraries
library(httr)  # Load data from Dropbox
require(spdep)
require(maptools)

library(ggplot2)
library(raster)
library(rgeos)

library(ggmap)

# install.packages("~/Downloads/rgdal_0.8-16.tgz")
library(rgdal)
library(rgeos)

# response <- GET(url="https://dl.dropboxusercontent.com/s/jl47y8w6da6x46m/Claudio.RData")
# load(rawConnection(response$content))
# rm(response)
load(paste0(filesDir,"/Claudio.Rdata"))
vlc.nb <- poly2nb(Carto,snap=0.01)

# setwd("~/docs/cursos/master_uv_bioestadistica/project/scripts/")
# load("./data/input_data.Rdata")

# Load shared functions and simulations
files2load <- c( "./utils.tfm.R", "./simu1.R", "./simu2.R", "./simu3.R")
for( file in files2load)
  source( file )

# Scenarios
simuList = c("simu1","simu2","simu3")

# Theta vector
v.theta = c(1.5,2,3)

# Scala factor to increase background expected cases
v.SF = c(1,2,4,10)

############
# Map
proj4string(Carto) <- CRS("+proj=utm +zone=30 +ellps=WGS84 +units=m +no_defs")
Carto@data$exp <- Eprostata

# res <- spTransform(Carto, CRS("+proj=longlat +datum=WGS84")) 
brks = quantile(Eprostata, seq(0, 1, length.out=256))
cols = colorRampPalette(c("#55FFFF", "grey10"))(255)
xlim = c(7.275e5,7.2751e5)
ylim=c(4.365e6,4.379e6)

####
# Plot with SPDEP package
plot(Carto,col=cols[findInterval(Eprostata, brks, all.inside=TRUE)],
     xlim=xlim,ylim=ylim)
title("Expected Risk values")

####
# Plot with SP package
xlim = c(720678,735386)
ylim = c(4362066,4383252)
spplot(Carto,"exp",xlim=xlim,ylim=ylim,main="Expected Risk values")

###########
# A little more elegantly with ggplot
Carto@data$exp <- Eprostata
xlim = c(720678,735386)
ylim = c(4362066,4383252)

# Create the clipping polygon
CP <- as(extent(xlim[1], xlim[2], ylim[1], ylim[2]), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(Carto))

# Clip the map
Carto.clip <- gIntersection(Carto, CP, byid=TRUE,id=rownames(Carto@data))

# Create the polygons
VLC.ggplot <- fortify(Carto.clip)
# Merge the Cases value to the polygons
VLC.ggplot <- cbind(VLC.ggplot, exp = Carto[VLC.ggplot$id, "exp"])

map.plot <- ggplot(data = VLC.ggplot) + geom_polygon(aes(x = long, y = lat, fill = exp, group = group)) + 
  coord_equal() + ggtitle("Expected values")  + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
map.plot

#  And now with a little cartographic thought.  As a choropleth with
# classes
VLC.ggplot$exp_cat <- cut(VLC.ggplot$exp, breaks = c(0, 1, 2, 5, 10, max(VLC.ggplot$exp)+1))
vlc.plot <- ggplot(data = VLC.ggplot) + geom_polygon(aes(x = long, y = lat, fill = exp_cat, 
                                            group = group)) + coord_equal() + scale_fill_brewer("Exp. val.", palette = "YlOrRd") +
 ggtitle("Expected values")  + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
vlc.plot

# And on google maps with the ggmap library
# VLC.ll <- spTransform(Carto, CRS("+proj=longlat +ellps=WGS84"))
# 
# vlc.map <- get_map(location = "Valencia, Spain", maptype = "roadmap", zoom = 12, color = "bw",source="osm")  #get map
# map.lim <- as.numeric(attr(vlc.map,"bb"))
# xlim = c(map.lim[2],map.lim[4])
# ylim = c(map.lim[1]-0.5,map.lim[3]+1)
# 
# CP <- as(extent(xlim[1], xlim[2], ylim[1], ylim[2]), "SpatialPolygons")
# proj4string(CP) <- CRS(proj4string(VLC.ll))
# 
# # Clip the map
# VLC.clip <- gIntersection(VLC.ll, CP, byid=TRUE,id=rownames(Carto@data))
# 
# # Create the polygons
# VLC.ggplot <- fortify(VLC.clip)
# # Merge the Cases value to the polygons
# VLC.ggplot <- cbind(VLC.ggplot,exp=Carto[VLC.ggplot$id, "exp"])
# VLC.ggplot$exp_cat <- cut(VLC.ggplot$exp, breaks = c(0, 1, 2, 5, 10, max(VLC.ggplot$exp)+1))
# 
# vlc.gmaps <- ggmap(vlc.map, extent = "panel") + geom_polygon(aes(x = long, y = lat, fill = exp_cat, group = group), data = VLC.ggplot, alpha = 0.7) + 
#   scale_fill_brewer("Exp. val.", palette = "YlOrRd")
# vlc.gmaps

##########
# Neighbours plot
######################################################## Create a
######################################################## neighborhood list
######################################################## from the
######################################################## shapefile
vlc_nb <- poly2nb(Carto.clip, queen = FALSE,snap=0.01)
summary(vlc_nb)
# Neighbour list object:
# Number of regions: 553 
# Number of nonzero links: 3040 
# Percentage nonzero weights: 0.9940845 
# Average number of links: 5.497288 
# Link number distribution:
#   
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  19 
# 1   5  36 130 155 108  62  18  22   4   4   1   1   2   2   1   1 
# 1 least connected region:
# 4625017006 with 1 link
# 1 most connected region:
# 4625015013 with 19 links

# Plot boundaries and edges
plot(Carto.clip, border = "grey60", axes = F)
title("Spatial Connectivity")
plot(vlc_nb, coordinates(Carto), pch = 19, cex = 0.6, add = TRUE)


############
# Scenario.1

obj <- simu1(Eprostata,5,1.5)
# plot(Carto,col=ifelse( Carto$CODIGO %in% obj$sel,"red4","navajowhite2" ),
#      xlim=xlim,ylim=ylim)
# title("Scenario.1 example")

VLC.ggplot$sel <- VLC.ggplot$id %in% obj$sel
vlc.plot <- ggplot(data = VLC.ggplot) + geom_polygon(aes(x = long, y = lat, fill = sel, group = group)) + 
  coord_equal() + scale_fill_brewer("Selected area", type="qual",palette=2) + 
  ggtitle("Scenario.1 example") + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
vlc.plot

# obj <- simu1(Eprostata,5,2,SF=2)
# plot(Carto,col=ifelse( Carto$CODIGO %in% obj$sel,"red4","navajowhite2" ),
#      xlim=xlim,ylim=ylim)
# title("Scenario.1 example")
# 
# obj <- simu1(Eprostata,5,3,SF=10)
# plot(Carto,col=ifelse( Carto$CODIGO %in% obj$sel,"red4","navajowhite2" ),
#      xlim=xlim,ylim=ylim)
# title("Scenario.1 example")

############
# Scenario.2
obj <- simu2(v.exp=Eprostata,all.nb=vlc.nb,theta=2,SF=1)
# plot(Carto,col=ifelse(Carto$CODIGO %in% obj$sel, "red3","navajowhite2" ),
#      xlim=xlim,ylim=ylim,)
# sel.idx <- unlist( lapply( unlist(obj$sel), FUN=function(x) which(names(obj$obs)==substr(x,6,10)) ) )
# title(main="Scenario.2 example",sub=paste0(round(sum(Eprostata[sel.idx])/sum(Eprostata)*100,2),"% expected cases over total cases in cluster."))

VLC.ggplot$sel <- VLC.ggplot$id %in% obj$sel
vlc.plot <- ggplot(data = VLC.ggplot) + geom_polygon(aes(x = long, y = lat, fill = sel, group = group)) + 
  coord_equal() + scale_fill_brewer("Selected area", type="qual",palette=2) + 
  ggtitle("Scenario.2 example") + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
vlc.plot

############
# Scenario.3
obj <- simu3( v.exp=Eprostata,
              all.nb=vlc.nb,
              theta=v.theta[3],NOC=20,SF=1,pct=0.01)
# plot(Carto,col=ifelse(Carto$CODIGO %in% unlist(obj$sel), "red3","navajowhite2" ),
#      xlim=xlim,ylim=ylim)
# title("Scenario.3 example")

VLC.ggplot$sel <- VLC.ggplot$id %in% unlist(obj$sel)
vlc.plot <- ggplot(data = VLC.ggplot) + geom_polygon(aes(x = long, y = lat, fill = sel, group = group)) + 
  coord_equal() + scale_fill_brewer("Selected area", type="qual",palette=2) + 
  ggtitle("Scenario.3 example") + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
vlc.plot

############
# BYM vs CROSSVAL
simu <- 3
theta <- 2
SF <- 2
rep <- 1
str.criteria <- "crossval"

plotResultMap <- function( simu, theta, SF, rep, str.criteria="crossval") {
    scen_results <- bivar.model[ bivar.model$simu==simu & bivar.model$theta==theta & bivar.model$SF==SF & bivar.model$rep==rep, ]
    
    scenFile <- paste0("./scenarios/simu",simu,"_",
                       gsub("([.])","",as.character(theta)),"_",as.character(SF),
                       "_r",as.character(rep),".dat")
    
    data <- dget( scenFile ) 
    sel.idx <- unlist( lapply( unlist(data$sel), FUN=function(x) which(names(data$obs)==substr(x,6,10)) ) )
    
    rownames(scen_results) <- names(data$obs)
    scen_results$code <- names(data$obs)
    
    # GGplot data preprocess
    Carto@data$exp <- Eprostata
    xlim = c(min(VLC.centers[rownames(VLC.centers)%in% unlist(data$sel), c("long")])- 1000, max(VLC.centers[rownames(VLC.centers) %in% unlist(data$sel),c("long")]) + 1000)
    ylim = c(min(VLC.centers[rownames(VLC.centers)%in% unlist(data$sel), c("lat")]) - 1000, max(VLC.centers[rownames(VLC.centers) %in% unlist(data$sel),c("lat")]) + 1000)
    
    # Create the clipping polygon
    CP <- as(extent(xlim[1], xlim[2], ylim[1], ylim[2]), "SpatialPolygons")
    proj4string(CP) <- CRS(proj4string(Carto))
    
    # Clip the map
    Carto.clip <- gIntersection(Carto, CP, byid=TRUE,id=as.character(Carto$CODIGO))
    
    # Create the polygons
    VLC.ggplot <- fortify(Carto.clip)
    # Merge the Cases value to the polygons
    VLC.ggplot <- cbind(VLC.ggplot, exp = Carto[VLC.ggplot$id, "exp"])
    
    VLC.ggplot$sel <- VLC.ggplot$id %in% unlist(data$sel)
    VLC.ggplot$criteria <- scen_results[ substr(VLC.ggplot$id,6,10) , str.criteria ]
    VLC.ggplot$criteria.sel <- VLC.ggplot$criteria * VLC.ggplot$sel 
      
    vlc.plot <- ggplot(data = VLC.ggplot) + geom_polygon(aes(x = long, y = lat, fill = criteria.sel, group = group)) + 
      coord_equal() + # scale_fill_brewer("Selected area", type="qual",palette=2) + 
      ggtitle(paste0("Scenario.",simu,";theta=",theta,";SF=",SF,";replica=",rep,"\n",str.criteria," result")) + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("") 
    
    # Add area labels   
    VLC.centers <- data.frame( coordinates(Carto) )
    VLC.centers$area <- substr(rownames(VLC.centers),6,10)
    colnames(VLC.centers)[1:2] <- c("long","lat")

 vlc.plot + geom_text(data=VLC.centers[rownames(VLC.centers)%in% unlist(data$sel),],aes(x=long,y=lat,label=area),size=3)

}

plotResultMap( simu, theta, SF, rep, str.criteria )
