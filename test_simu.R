# test_simu.R
################################
# Date: 10/02/2014
# Author: CNM 
# Project: TFM
# Description: script to test simu1, simu2, simu3 functions


# Clean workspace & set wd
rm(list=ls())
setwd("~/GitHub/tfm")

# Libraries
library(httr)  # Load data from Dropbox
require(spdep)
require(maptools)

response <- GET(url="https://dl.dropboxusercontent.com/s/jl47y8w6da6x46m/Claudio.RData")
load(rawConnection(response$content))
rm(response)
# setwd("~/docs/cursos/master_uv_bioestadistica/project/scripts/")
# load("./data/input_data.Rdata")

# Load shared functions and simulations
files2load <- c( "./utils.tfm.R", "./simu1.R", "./simu2.R", "./simu3.R")
for( file in files2load)
  source( file )

###################
# Testing

# Definicion del vector de thetas
v.theta = c(1.5,2,3)

# Factor de escala para la los valores esperados
v.SF = c(1,2,4,10)

# nb object from spatial polygons.
vlc.nb <- poly2nb(Carto,snap=0.01)


############
# Scenario.1
brks = quantile(Eprostata, seq(0, 1, length.out=256))
cols = colorRampPalette(c("#55FFFF", "grey10"))(255)
xlim = c(7.275e5,7.2751e5)
ylim=c(4.365e6,4.379e6)

plot(Carto,col=cols[findInterval(Eprostata, brks, all.inside=TRUE)],
     xlim=xlim,ylim=ylim)

obj <- simu1(Eprostata,5,1.5)
plot(Carto,col=ifelse( Carto$CODIGO %in% obj$sel,"red4","navajowhite2" ),
     xlim=xlim,ylim=ylim)

obj <- simu1(Eprostata,5,2,SF=2)
plot(Carto,col=ifelse( Carto$CODIGO %in% obj$sel,"red4","navajowhite2" ),
     xlim=xlim,ylim=ylim)

obj <- simu1(Eprostata,5,3,SF=10)
plot(Carto,col=ifelse( Carto$CODIGO %in% obj$sel,"red4","navajowhite2" ),
     xlim=xlim,ylim=ylim)


############
# Scenario.2
obj <- simu2(v.exp=Eprostata,all.nb=vlc.nb,theta=2,SF=1)
plot(Carto,col=ifelse(Carto$CODIGO %in% obj$sel, "red3","navajowhite2" ),
     xlim=xlim,ylim=ylim)


############
# Scenario.3
obj <- simu3( v.exp=Eprostata,
              all.nb=vlc.nb,
              theta=v.theta[3],NOC=20,SF=1,pct=0.01)
# Plot
plot(Carto,col=ifelse(Carto$CODIGO %in% unlist(obj$sel), "red3","navajowhite2" ),
     xlim=xlim,ylim=ylim)
