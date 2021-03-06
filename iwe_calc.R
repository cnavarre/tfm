# iwe_cacl.R
################################
# Date: 29/03/2014
# Author: CNM 
# Project: TFM
# Description: Importance weigthed estimates

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

require(spdep)
require(maptools)

# Load shared functions and simulations
files2load <- c( "./utils.tfm.R", "./simu1.R", "./simu2.R", "./simu3.R")
for( file in files2load)
  source( file )

# Load map data
# response <- GET(url="https://www.dropbox.com/s/xtbhm2czs5ox5es/Claudio.Rdata")
# load(rawConnection(response$content))
# rm(response)
load(paste0(filesDir,"/Claudio.Rdata"))
vlc.nb <- poly2nb(Carto,snap=0.01)

###################
# Constants
# Number of regions
Q = length(vlc.nb)

# Definicion del vector de thetas
v.theta = c(1.5,2,3)

# Factor de escala para la los valores esperados
v.SF = c(1,2,4,10)

# Threshold for ROC classification
threshold <- 0.8

# DEBUG START
# simu <- "simu1"
# theta = 1.5
# SF = 10
# theta.str = gsub("[.]","",as.character(theta))
# DEBUG END
simuList = c("simu1","simu2","simu3")

iwe.table <- data.frame()

# Result file Path
resultsPath <- paste0(filesDir,"/output/")

for( simu in simuList ) {
for( theta in v.theta) {
for( SF in v.SF) {
    # Print new parameters
    print(paste("> Simu=",simu),quote=F)
    print(paste("> Theta=",theta),quote=F)
    print(paste("> SF=",SF),quote=F)
  
  
    theta.str = gsub("[.]","",as.character(theta))
    
    pattern.str <- paste0(simu,"_",theta.str,"_", SF,"_")
    resultsList <- paste0(resultsPath,dir(path=resultsPath,pattern=pattern.str))    
    
    for( resultFile in resultsList) {
      # DEBUG
      # resultFile <- resultsList[4]
      # DEBUG END
      
      # Print filename
      filename <- basename(resultFile)
      print(paste(">> Filename: ", filename), quote=F)

      bym <- readRDS(file=resultFile)
      bym.mcmc <- as.mcmc(bym)
      
      scenFile <- paste0("./scenarios/",sub("\\.[[:alnum:]]+$", "", basename(resultFile)),".dat")
      data <- dget( scenFile )
      sel.idx <- unlist( lapply( unlist(data$sel), FUN=function(x) which(names(data$obs)==substr(x,6,10)) ) )
      O <- data$obs # Y_i

      # Calculate E_i via original value + SF
      E <- Eprostata * SF
      
      # Using same theta.post calculated through BYM model
#       theta.post <- t(t(bym.mcmc$sims.matrix[,grep("R\\[",colnames(bym.mcmc$sims.matrix))]) * E )
      # Number of replicates
      K <- bym.mcmc$n.sims # dim(theta.post)[1]

      # Calculating an approximate sample of theta.rep
      m <- bym.mcmc$sims.matrix[,grep("m",colnames(bym.mcmc$sims.matrix))]
      sp <- bym.mcmc$sims.matrix[,grep("sp\\[",colnames(bym.mcmc$sims.matrix))]
      sdsp <- bym.mcmc$sims.matrix[,grep("sdsp",colnames(bym.mcmc$sims.matrix))]
      sdhet <- bym.mcmc$sims.matrix[,grep("sdhet",colnames(bym.mcmc$sims.matrix))]
      
      # sp.rep_ij ~ N( sum_k_in_vlc.nb[[i]]( sp[k,j] ) / n_i , sdsp[j] / n_i ) with n_i=number of neigbours of region i
      sp.rep <- sapply(1:Q,function(i) sapply(1:K, function(j) rnorm( 1, mean=mean(sp[j,vlc.nb[[i]] ]) , sd=sdsp[j] / card(vlc.nb)[i] ) ) )
      het.rep <- sapply( 1:K, function(j) rnorm(1,mean=0,sd=sdhet[j]) )
      
      # theta.rep_i = exp(m_i + het.rep_i + sp.rep_i )
      theta.rep <- exp( m + sp.rep + het.rep )
      theta.post <- t(t(theta.rep) * E)
            
      # w\i(\theta^[k])
      w <- 1 / matrix(apply(theta.post,1,function(x){dpois(O,x)}),nrow=nrow(theta.post),byrow=T)
      # P.rep      
      P.rep <- t( apply(theta.post,1, function(theta.row) { ppois(O, theta.row) - 0.5 * dpois(O,theta.row) } ) )
      # Prob(Y.rep<=O_i|O\_i)
      P.Y.rep <- apply(P.rep * w,2,sum) / apply(w,2,sum)
      
      # Read P.Y.rep
      P.Y.rep <- dget( file=paste0(filesDir, "/results/crossval/","pyrep_",basename(resultFile) ) )

      # Calculate SPEC. y SENS.
      iwe.table <- rbind(iwe.table,data.frame( file=filename, simu = substr(simu,5,5),
                                               rep= as.numeric( substr( basename(resultFile), regexpr("r",filename)+1, regexpr("[.]",filename)-1 ) ),
                                               theta, SF,
                                               t(roc.curve( obs = 1:Q %in% sel.idx, pred=P.Y.rep,th=threshold )),check.names=F))
 
      # Write P.Y.rep into a result file
      # dput( x=P.Y.rep, file=paste0(filesDir, "/results/crossval/","pyrep_",basename(resultFile) ) )
    } # END for( resultFile in resultsList) 
    
} # end for( SF in v.SF)
} # end for( theta in v.theta)
} # end for( simu in simuList ) 

# Save table
# iwe.table$rep <- as.numeric(with(iwe.table,substr(file,regexpr("r",file)+1, regexpr("[.]",file)-1) ))
save(iwe.table,file="./iwe_table.RData")
load("iwe_table.RData")
# Plots
# colnames(iwe.table)[4:5] <- c("FPR","TPR")
ggplot(iwe.table,aes(y=1-FPR,x=as.factor(theta),fill=as.factor(theta))) + geom_boxplot() + facet_grid(simu~SF) + ggtitle("IWE method \nspec.(1-FPR)")
ggplot(iwe.table,aes(y=TPR,x=as.factor(theta),fill=as.factor(theta))) + geom_boxplot() + facet_grid(simu~SF) + ggtitle("IWE method \nTPR (sens.)")

ggplot(iwe.table,aes(y=TPR,x=1-FPR,colour=as.factor(theta),xmin=0,xmax=1,ymin=0,ymax=1)) + geom_point() + facet_grid(simu~SF) + ggtitle("IWE method \nSensitivity vs Specificity")


# Tables
# dcast(iwe.table,simu~theta+SF,function(x) mean(x,na.rm=T),value.var="TPR")
# dcast(iwe.table,simu~theta+SF,function(x) mean(x,na.rm=T),value.var="FPR")
