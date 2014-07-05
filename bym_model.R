# bym_model.R
################################
# Date: 10/02/2014
# Author: CNM 
# Project: TFM
# Description: 

# Clean workspace & set wd
rm(list=ls())
setwd("~/GitHub/tfm")

# Libraries
library(httr)  # Load data from Dropbox
require(spdep)
require(maptools)
library(R2WinBUGS)

response <- GET(url="https://dl.dropboxusercontent.com/s/jl47y8w6da6x46m/Claudio.RData")
load(rawConnection(response$content))
rm(response)

# Prepare winbugs data
carto.nb <- poly2nb(Carto,snap=1)
carto.wb <- nb2WB(carto.nb)
rm(list=c("Carto","Eprostata","carto.nb"))

# Load scenario
scenList <- paste0("./scenarios/",dir(path="scenarios",pattern="simu1"))
# scenList <- paste0("./scenarios/",dir(path="scenarios",pattern="simu2"))
# scenList <- paste0("./scenarios/",dir(path="scenarios",pattern="simu3"))

iters = 3500
burn = 500

for( scen in scenList ) {
  # path <- dirname( strScenario )
  # filename <- basename( strScenario )
  outputfile <- paste0("output","/",sub("\\.[[:alnum:]]+$", "", basename(scen)),".Rda")
  data <- dget( scen )

  O <- data$obs
  E <- data$exp.base # Expected values modified

  # Modelo Suavizado Besag York y Mollie
  bym.data <- list( O=O, E=E, n=length(O), 
                    adj = carto.wb$adj, weights=carto.wb$weights,
                    num = carto.wb$num )
  bym.inits <- function() {list(prechet = 1, precsp = 1,m = 0,
                                het = rep(0,length(O)), sp = rep(0,length(O)))}
  bym.params <- c("m","sdhet","sdsp","R","sp")
  model.file <- paste(getwd(),"/bym.model.bugs",sep="")
  bym <- bugs( data = bym.data, inits = bym.inits, parameters=bym.params,
               model.file = model.file, n.chains = 3, 
               working.directory="~/.wine/drive_c/temp/Rtmp/", clearWD=T,
               bugs.directory="c:/Program Files (x86)/WinBUGS14/",
               n.iter = iters, n.burnin = burn, n.thin = 1,
               debug = F, DIC = F)

  saveRDS( object=bym,file=outputfile )

  # bym <- readRDS(file=outputfile)

  # bym.mcmc <- as.mcmc(bym)
  # bym.m <- as.matrix(bym.mcmc)
}