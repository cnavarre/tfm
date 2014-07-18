# sens_post_table.R
################################
# Author: CNM 
# Date: 13/07/2014
# Project: TFM
# Description: Posterior analysis of BYM model in cluster detection. 
#              SENSITIVITY.

# Clean workspace & set wd
rm(list=ls())
# filesDir <- "C:/Users/claudio/Dropbox/Claudio-Migue"
# workDir <- "~/GitHub/tfm"
# filesDir <- "/host/Users/claudio/Dropbox/Claudio-Migue"
# workDir <- "/host/Users/claudio/Documents/GitHub/tfm"
filesDir <- "~/Dropbox/Claudio-Migue"
workDir <- "~/Documents/docs/projects/tfm"
setwd(workDir)

# Libraries
library(httr)  # Load data from Dropbox
library(R2WinBUGS)
library(coda)
library(xtable)
library(ggplot2)
library(reshape2)
library(spdep)

# Load shared functions and simulations
files2load <- c( "./utils.tfm.R", "./simu1.R", "./simu2.R", "./simu3.R")
for( file in files2load)
  source( file )

# Load map data
response <- GET(url="https://dl.dropboxusercontent.com/s/jl47y8w6da6x46m/Claudio.RData")
load(rawConnection(response$content))
rm(response)
vlc.nb <- poly2nb(Carto,snap=0.01)


###################
# Constants
# Number of regions
Q = length(vlc.nb)

# Definicion del vector de thetas
v.theta = c(1.5,2,3)

# Factor de escala para la los valores esperados
v.SF = c(1,2,4,10)

# DEBUG START
theta = 1.5
SF = 10
theta.str = gsub("[.]","",as.character(theta))
# DEBUG END
simuList = c("simu1","simu2","simu3")

# Result file Path
# resultsPath <- paste0(filesDir,"/output/")  # old data without sp draws
resultsPath <- paste0(filesDir,"/output2/") 

# Result table with FPR (1-spec) + PDR (sens) for model BYM based on D(0.8,1)
bym.table <- data.frame()

for( simu in simuList ) {
  for( theta in v.theta) {
    for( SF in v.SF) {
      theta.str = gsub("[.]","",as.character(theta))
      
      pattern.str <- paste0(simu,"_",theta.str,"_", SF,"_")
      resultsList <- paste0(resultsPath,dir(path=resultsPath,pattern=pattern.str))    
            
      for( resultFile in resultsList) {
        bym <- readRDS(file=resultFile)      
        bym.mcmc <- as.mcmc(bym)
        bym.mcmc$sims.matrix <- bym.mcmc$sims.matrix[,-(1:3)] # We are interested only in posterior dist. of \theta_i
        
        scenFile <- paste0("./scenarios/",sub("\\.[[:alnum:]]+$", "", basename(resultFile)),".dat")
        data <- dget( scenFile )
        #       data$SMR <- data$obs / data$exp
        
        sel.idx <- unlist( lapply( unlist(data$sel), FUN=function(x) which(names(data$obs)==substr(x,6,10)) ) )
        
        # Classify with decission rule function
        filename <- basename(resultFile)
        bym.table <- rbind(bym.table,data.frame( file = filename, simu= as.numeric(substr(simu,5,5)),
                                                 rep = as.numeric( substr( basename(resultFile), regexpr("r",filename)+1, regexpr("[.]",filename)-1 ) ),
                                                 theta, SF, 
                                                 t( roc.curve( seq(Q) %in% sel.idx, 
                                                    destring(decision.rule(bym.mcmc$sims.matrix[,col_name("R",seq(Q))] )),0.5 ) ) , check.names=F ) )
      } # END for( resultFile in resultsList )
        
    } # for( SF in v.SF) 
  } # END for( theta in v.theta) 
} # END for( simu in simuList ) 

# Save table
save(bym.table,file="./spec_sens_table.RData")

# Plots
ggplot(bym.table,aes(y=FPR,x=as.factor(theta),fill=as.factor(theta))) + geom_boxplot() + facet_grid(simu~SF) + ggtitle("BYM method \nFPR (1-spec.)")
ggplot(bym.table,aes(y=TPR,x=as.factor(theta),fill=as.factor(theta))) + geom_boxplot() + facet_grid(simu~SF) + ggtitle("BYM method \nTPR (sens.)")

ggplot(bym.table,aes(y=TPR,x=FPR,colour=as.factor(theta),xmin=0,xmax=1,ymin=0,ymax=1)) + geom_point() + facet_grid(~SF) + ggtitle("BYM method \nTPR vs FPR")
