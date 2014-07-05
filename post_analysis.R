# post_analysis.R
################################
# Date: 24/02/2014
# Author: CNM 
# Project: TFM
# Description: Posterior analysis of BYM model in cluster detection

# Clean workspace & set wd
rm(list=ls())
# filesDir <- "C:/Users/claudio/Dropbox/Claudio-Migue"
# workDir <- "~/GitHub/tfm"
filesDir <- "/host/Users/claudio/Dropbox/Claudio-Migue"
workDir <- "/host/Users/claudio/Documents/GitHub/tfm"
setwd(workDir)

# Libraries
library(httr)  # Load data from Dropbox
library(R2WinBUGS)
library(coda)
library(xtable)
library(ggplot2)
library(reshape2)

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

# DEBUG START
theta = 1.5
SF = 10
theta.str = gsub("[.]","",as.character(theta))
# DEBUG END
simuList = c("simu1","simu2","simu3")

# Result file Path
# resultsPath <- paste0(filesDir,"/output/")  # old data without sp draws
resultsPath <- paste0(filesDir,"/output2/") 

rem(c("table1","table2","table3","table4","table5","table6","table7"))

table7 <- data.frame(false.positive = numeric(), 
                     simu = character(), 
                     theta = numeric(), SF=numeric(),stringsAsFactors=FALSE)
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
      {  }
      )
      
      
      # Sensitivity 
      # Classifying cluster areas
#       switch( simu,
#        simu1={ifelse( exists("class.clust"),
#                 class.clust <- rbind(class.clust,destring(decision.rule( bym.mcmc$sims.matrix[ ,sel.idx + 3 ]) )),
#                 class.clust <- destring(decision.rule( bym.mcmc$sims.matrix[ ,sel.idx + 3 ]) )
#               )},
#        {
#          ifelse( exists("class.clust"),
#                  class.clust <- c(class.clust,mean(destring(decision.rule( bym.mcmc$sims.matrix[ ,sel.idx + 3 ]) ) ) ),
#                  class.clust <- mean(destring(decision.rule( bym.mcmc$sims.matrix[ ,sel.idx + 3 ]) ) ) )
#        } )
       
      # False positive rate (1-specificity) 
      # Classifying background areas
      ifelse( exists("class.bg"),
              class.bg <- c(class.bg,mean(destring(decision.rule( bym.mcmc$sims.matrix[ ,-(sel.idx + 3) ]) ) ) ),
              class.bg <- mean(destring(decision.rule( bym.mcmc$sims.matrix[ ,-(sel.idx + 3) ]) ) ) ) 
      
#       # SMR of clusters areas
#       switch( simu,
#       simu1={ ifelse( exists("post.mean.clust"),
#                       post.mean.clust <- rbind(post.mean.clust, apply(bym.mcmc$sims.matrix[ ,sel.idx + 3 ],2,mean) ), 
#                       post.mean.clust <- apply(bym.mcmc$sims.matrix[ ,sel.idx + 3 ],2,mean) )              
#       },
#       simu2= {
#         ifelse( exists("post.mean.clust"),
#                 post.mean.clust <- rbind(post.mean.clust, mean(bym.mcmc$sims.matrix[ ,sel.idx + 3 ]) ), 
#                 post.mean.clust <- mean(bym.mcmc$sims.matrix[ ,sel.idx + 3 ]) )
#       },
#       simu3= {
#         ifelse( exists("post.mean.clust"),
#                 post.mean.clust <- rbind(post.mean.clust, mean(bym.mcmc$sims.matrix[ ,sel.idx + 3 ]) ), 
#                 post.mean.clust <- mean(post.mean.clust$sims.matrix[ ,sel.idx + 3 ]) )
#       } )
      
#       # SMR of background areas
#       ifelse( exists("bym.m.bg"),
#               bym.m.bg <- rbind(bym.m.bg, bym.mcmc$sims.matrix[ ,-(sel.idx + 3) ] ), 
#               bym.m.bg <- bym.mcmc$sims.matrix[ ,-(sel.idx + 3) ] )
    }
     
#     ###########
#     # Table.1. Posterior mean relative risk estimates for "raised-risk" areas.
#     switch( simu,
#             simu1={
#                   colnames( post.mean.clust ) <- c("0-20% area", "20-40% area", "40-60% area", "60-80% area", "80-100% area")
#                   
#                   
#                   colnames( table1 ) <- c( colnames( table1 )[-dim(table1)[2]], paste0("\theta=",theta,";SF=",SF) )
#             },
#             simu2={   
#               ifelse( exists("table2"),
#                       table2 <- cbind(table2, mean(post.mean.clust) ),
#                       table2 <- data.frame( mean(post.mean.clust) ) )    
#             
#               colnames( table2 ) <- c( colnames( table2 )[-dim(table2)[2]], paste0("\theta=",theta,";SF=",SF) )
#               rownames( table2 ) <- "1% cluster"
#             },
#             simu3={
#               ifelse( exists("table3"),
#                       table3 <- cbind(table3, mean(post.mean.clust) ),
#                       table3 <- data.frame( mean(post.mean.clust) ) )    
#               
#               colnames( table3 ) <- c( colnames( table3 )[-dim(table3)[2]], paste0("\theta=",theta,";SF=",SF) )
#               rownames( table3 ) <- "20x1% clusters"
#             }
#     )

    ###########
    # Table.[4,5,6]. Sensitivity: 1-False.negative rate
#     switch( simu,
#         simu1={
#           ifelse( exists("table4"),
#                   table4 <- cbind(table4, apply(class.clust,2,mean) ),
#                   table4 <- data.frame( sensitivity=apply(class.clust,2,mean) ) )    
#           
#           rownames( table4 ) <- c("0-20% area", "20-40% area", "40-60% area", "60-80% area", "80-100% area")
#           colnames( table4 ) <- c( colnames( table4 )[-dim(table4)[2]], paste0("\theta=",theta,";SF=",SF) )
#         },
#         simu2={   
#           ifelse( exists("table5"),
#                   table5 <- cbind(table5, mean(class.clust) ),
#                   table5 <- data.frame( mean(class.clust )) )    
#           
#           colnames( table5 ) <- c( colnames( table5 )[-dim(table5)[2]], paste0("\theta=",theta,";SF=",SF) )
#           rownames( table5 ) <- "1% cluster"
#         },
#         simu3={
#           ifelse( exists("table6"),
#                   table6 <- cbind(table6, mean(class.clust )),
#                   table6 <- data.frame( mean(class.clust)) )    
#           
#           colnames( table6 ) <- c( colnames( table6 )[-dim(table6)[2]], paste0("\theta=",theta,";SF=",SF) )
#           rownames( table6 ) <- "20x1% clusters"
#         }
#       )
    
    ###########
    # Table.7. False positive rates: 1-specificity    
    table7 <- rbind(table7, data.frame(false.positive=mean(class.bg),simu,theta,SF) )
    
    ###########
    # Figure.1. Histograms raw SMRs
#     ifelse( exists("fig1.df"),
#       fig1.df <- rbind(fig1.df,data.frame(SMR=data$SMR[-sel.idx])),
#       fig1.df <- data.frame(SMR=data$SMR[-sel.idx]) )
#       
     # rm(list=c("bym.m.clust")) #,bym.m.bg))
}
}
}

# print(xtable(table1),type="html")
# save(table1,table2,table3,file="tables.Rdata")
load("tables.Rdata")
table1
table2
table3

# Sensitivity (Prob. of true positive)
# save(table4,table5,table6,file="tables2.Rdata")
load("tables2.Rdata")
sum.table4 <- apply(table4,2,mean)
sum.table4
table5
table6

# False positive rates (1-specificity)
# save(table7,file="tables3.Rdata")
load("tables3.Rdata")
table7.b <- dcast(table7,simu~theta+SF,value.var="false.postive")

# Expected mean for Simu1 0-20%,20-40%,40-60%,60-80%,80-100% areas
hist(Eprostata,breaks=60,xlab="Expected cases",main="Histogram")
points(cbind(expval,0),pch=4,col="red")
