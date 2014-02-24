# post_analysis.R
################################
# Date: 24/02/2014
# Author: CNM 
# Project: TFM
# Description: Posterior analysis of BYM model in cluster detection

# Clean workspace & set wd
rm(list=ls())
filesDir <- "C:/Users/claudio/Dropbox/Claudio-Migue"
workDir <- "~/GitHub/tfm"
setwd(workDir)

# Libraries
library(httr)  # Load data from Dropbox
library(R2WinBUGS)
library(coda)
library(xtable)

###################
# Constants
# Definicion del vector de thetas
v.theta = c(1.5,2,3)

# Factor de escala para la los valores esperados
v.SF = c(1,2,4,10)

# DEBUG START
theta = 3
SF = 1
theta.str = gsub("[.]","",as.character(theta))

simu = "simu1"
# DEBUG END


# Result file Path
resultsPath <- paste0(filesDir,"/output/")


for( theta in v.theta) {
  for( SF in v.SF) {
    theta.str = gsub("[.]","",as.character(theta))
    
    pattern.str <- paste0(simu,"_",theta.str,"_", SF,"_")
    resultsList <- paste0(resultsPath,dir(path=resultsPath,pattern=pattern.str))    
    
    rm(bym.m)
    for( resultFile in resultsList) {
      bym <- readRDS(file=resultFile)      
      bym.mcmc <- as.mcmc(bym)
      
      scenFile <- paste0("./scenarios/",sub("\\.[[:alnum:]]+$", "", basename(resultFile)),".dat")
      data <- dget( scenFile )
      
      sel.idx <- unlist( lapply( 1:5, FUN=function(x) which(names(data$obs)==substr(data$sel[x],6,10)) ) )
      
      ifelse( exists("bym.m"),
          bym.m <- rbind(bym.m, bym.mcmc$sims.matrix[ ,sel.idx + 3 ] ), 
          bym.m <- bym.mcmc$sims.matrix[ ,sel.idx + 3 ] )
    }
     
    colnames( bym.m ) <- c("0-20% area", "20-40% area", "40-60% area", "60-80% area", "80-100% area")

    ifelse( exists("table1"),
            table1 <- cbind(table1, apply(bym.m,2,mean) ),
            table1 <- data.frame( apply(bym.m,2,mean) ) )    
    
    colnames( table1 ) <- c( colnames( table1 )[-dim(table1)[2]], paste0("\theta=",theta,";SF=",SF) )
    
    rm(bym.m)
  }
}
