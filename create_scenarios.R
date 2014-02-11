# create_scenarios.R
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
require("maptools")

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
# Constants
# Definicion del vector de thetas
v.theta = c(1.5,2,3)

# Factor de escala para la los valores esperados
v.SF = c(1,2,4,10)

# nb object from spatial polygons.
vlc.nb <- poly2nb(Carto,snap=1)

###################
# Create scenarios
for( theta in v.theta)
  for( SF in v.SF) {
    # Simu1
    obj <- simu1(Eprostata,5,theta=theta,SF=SF)
    dput(x=obj,
         file=paste0("./scenarios/simu1_",
                     gsub("([.])","",as.character(theta)),"_",as.character(SF),".dat") )
    
    # Simu2
    obj <- simu2(Eprostata,5,theta=theta,SF=SF,all.nb=vlc.nb)
    dput(x=obj,
         file=paste0("./scenarios/simu2_",
                     gsub("([.])","",as.character(theta)),"_",as.character(SF),".dat") )
    
    # Simu3
    obj <- simu3(Eprostata,5,theta=theta,SF=SF,all.nb=vlc.nb)
    dput(x=obj,
         file=paste0("./scenarios/simu3_",
                     gsub("([.])","",as.character(theta)),"_",as.character(SF),".dat") )
  }
   