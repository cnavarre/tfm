# smr_post_table.R
################################
# Author: CNM 
# Date: 13/07/2014
# Project: TFM
# Description: Posterior analysis of BYM model in cluster detection. 
#              TABLE of SMR for increased risk areas.

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
# Definicion del vector de thetas
v.theta = c(1.5,2,3)

# Factor de escala para la los valores esperados
v.SF = c(1,2,4,10)

# DEBUG START
# simu = "simu1"
# theta = 1.5
# SF = 10
# DEBUG END
simuList = c("simu1","simu2","simu3")

# Result file Path
# resultsPath <- paste0(filesDir,"/output/")  # old data without sp draws
resultsPath <- paste0(filesDir,"/output2/") 

rem(c("table1","table2","table3"))

Sys.time()
for( simu in simuList ) {
for( theta in v.theta) {
for( SF in v.SF) {
      theta.str = gsub("[.]","",as.character(theta))
      
      pattern.str <- paste0(simu,"_",theta.str,"_", SF,"_")
      resultsList <- paste0(resultsPath,dir(path=resultsPath,pattern=pattern.str))    
      
      rem(c("post.mean.clust","class.clust","class.bg"))
      
      for( resultFile in resultsList) {
        bym <- readRDS(file=resultFile)      
        bym.mcmc <- as.mcmc(bym)
        bym.mcmc$sims.matrix <- bym.mcmc$sims.matrix[,-(1:3)] # We are interested only in posterior dist. of \theta_i
        
        scenFile <- paste0("./scenarios/",sub("\\.[[:alnum:]]+$", "", basename(resultFile)),".dat")
        data <- dget( scenFile )
        #       data$SMR <- data$obs / data$exp
        
        sel.idx <- unlist( lapply( unlist(data$sel), FUN=function(x) which(names(data$obs)==substr(x,6,10)) ) )
        switch( simu,
                simu1={ ifelse( exists("expval"),
                                expval <- apply(rbind(expval,Eprostata[sel.idx]),2,mean),
                                expval <- Eprostata[sel.idx])
                        names(expval) <- c("0-20% area", "20-40% area", "40-60% area", "60-80% area", "80-100% area")
                },
            {  } )
        
        # SMR of clusters areas
        switch( simu,
        simu1={ ifelse( exists("post.mean.clust"),
                        post.mean.clust <- rbind(post.mean.clust, apply(bym.mcmc$sims.matrix[ ,col_name( "R",sel.idx) ],2,mean) ), 
                        post.mean.clust <- apply(bym.mcmc$sims.matrix[ , col_name("R", sel.idx ) ],2,mean) )              
        },
        simu2= {
          ifelse( exists("post.mean.clust"),
                  post.mean.clust <- rbind(post.mean.clust, mean(bym.mcmc$sims.matrix[ , col_name( "R", sel.idx ) ]) ), 
                  post.mean.clust <- mean(bym.mcmc$sims.matrix[ , col_name( "R", sel.idx ) ]) )
        },
        simu3= {
          ifelse( exists("post.mean.clust"),
                  post.mean.clust <- rbind(post.mean.clust, mean(bym.mcmc$sims.matrix[ , col_name( "R", sel.idx ) ]) ), 
                  post.mean.clust <- mean(bym.mcmc$sims.matrix[ , col_name( "R", sel.idx ) ]) )
        } )
  
        # SMR of background areas
#         ifelse( exists("bym.m.bg"),
#                 bym.m.bg <- rbind(bym.m.bg, bym.mcmc$sims.matrix[ ,-(sel.idx + 3) ] ), 
#                 bym.m.bg <- bym.mcmc$sims.matrix[ ,-(sel.idx + 3) ] )
      } # for( resultFile in resultsList )
      
        
        
#     ###########
#     # Table.1. Posterior mean relative risk estimates for "raised-risk" areas.
      switch( simu,
            simu1={
                tmp.df <- data.frame( post.mean.clust, theta=theta, SF=SF,row.names=NULL )
                colnames(tmp.df) <- c("0-20% area", "20-40% area", "40-60% area", "60-80% area", "80-100% area","theta","SF") ;
              
                ifelse( exists("table1") ,
                      table1 <- rbind(table1, tmp.df ) ,                      
                      table1 <- tmp.df )
                  
#                   colnames( table1 ) <- c( colnames( table1 )[-dim(table1)[2]], paste0("\theta=",theta,";SF=",SF) )
            },
            simu2={   
              ifelse( exists("table2"),
                      table2 <- cbind(table2, mean(post.mean.clust) ),
                      table2 <- data.frame( mean(post.mean.clust) ) )    
            
              colnames( table2 ) <- c( colnames( table2 )[-dim(table2)[2]], paste0("\theta=",theta,";SF=",SF) )
              rownames( table2 ) <- "1% cluster"
            },
            simu3={
              ifelse( exists("table3"),
                      table3 <- cbind(table3, mean(post.mean.clust) ),
                      table3 <- data.frame( mean(post.mean.clust) ) )    
              
              colnames( table3 ) <- c( colnames( table3 )[-dim(table3)[2]], paste0("\theta=",theta,";SF=",SF) )
              rownames( table3 ) <- "20x1% clusters"
            }
          )
} # end for( SF in v.SF)
} # end for( theta in v.theta) 
} # end for( simu in simuList )
Sys.time()

# print(xtable(table1),type="html")
save(table1,table2,table3,file="tables.Rdata")
load("tables.Rdata")
table1
table2
table3

# TABLA.1: 
t(aggregate(table1[,-7:-6],by=list(SF=table1$SF,theta=table1$theta),mean))
