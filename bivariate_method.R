# bivariate_method.R
################################
# Date: 29/03/2014
# Author: CNM 
# Project: TFM
# Description: Different criteria selection method to classify clusters areas

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

bivar.model <- data.frame()

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
      # resultFile <- resultsList[1]
      # DEBUG END
      
      # Print filename
      filename <- basename(resultFile)
      print(paste(">> Filename: ", filename),quote=F )
      
      # Get replicate number
      rep <- as.numeric( substr( filename, regexpr("r",filename)+1, regexpr("[.]",filename)-1 ) )
      
      # Read P.Y.rep: P(Y_rep <= O_i| O_{j\i})
      P.Y.rep <- dget( file=paste0(filesDir, "/results/crossval/","pyrep_v2_",basename(resultFile) ) )

      # Read theta_post (BYM)
      bym <- readRDS(file=resultFile)
      bym.mcmc <- as.mcmc(bym)
      theta_post <- bym.mcmc$sims.matrix[,col_name("R",seq(Q))] # Posterior dist. of \theta
      
      # Get increased areas indexes
      scenFile <- paste0("./scenarios/",sub("\\.[[:alnum:]]+$", "", basename(resultFile)),".dat")
      data <- dget( scenFile ) 
      sel.idx <- unlist( lapply( unlist(data$sel), FUN=function(x) which(names(data$obs)==substr(x,6,10)) ) )
      
      # Calculate ( X_i, Y_i ) with X_i = P( theta_post_i > 1 ) and Y_i = P( Y_rep <= O_i|O_{j\i} )
      bivar.model <- rbind( bivar.model, data.frame( bym=apply(theta_post,2,function(x) sum(x>1)/length(x)), 
                                                     crossval=P.Y.rep, 
                                                     sel=as.numeric(seq(Q) %in% sel.idx),
                                                     order=match(seq(Q),sel.idx,nomatch=0),
                                                     theta=theta,
                                                     rep=rep,
                                                     simu=as.numeric(substr(simu,5,5)),
                                                     SF=SF) )
    }
    
} # end for( SF in v.SF)
} # end for( theta in v.theta)
} # end for( simu in simuList ) 

rownames( bivar.model ) <- seq(nrow(bivar.model))

head(bivar.model)
dim(bivar.model)

# Calculate mean ,max and wmean selection criteria 
bivar.model$mean <- apply( bivar.model[,1:2], 1, mean, na.rm=T )
bivar.model$max <- apply( bivar.model[,1:2], 1, max, na.rm=T )
bivar.model$wmean <- 0


for( simu in simuList) {
for( theta in v.theta) {
for( SF in v.SF) {
    print(paste("> Simu=",simu),quote=F)
    print(paste("> Theta=",theta),quote=F)
    print(paste("> SF=",SF),quote=F)
for( rep in 1:10) {
    print(paste("> rep=",rep),quote=F)
    temp <- bivar.model[ bivar.model$simu==as.numeric( substr(simu,5,5) ) & 
                           bivar.model$theta==theta & 
                           bivar.model$SF==SF & 
                           bivar.model$rep==rep,]
    # Calculamos los pesos
    ww <- svd(var(temp[,1:2],na.rm=T))$u[,1]
    
    if( prod(ww) > 0 ){
      bivar.model[ bivar.model$simu==as.numeric( substr(simu,5,5) ) & 
                     bivar.model$theta==theta & 
                     bivar.model$SF==SF & 
                     bivar.model$rep==rep,]$wmean <- apply( temp[,c("bym","crossval")], 1, function(v) weighted.mean(v,ww) )
    }else{
      bivar.model[ bivar.model$simu==as.numeric( substr(simu,5,5) ) & 
                     bivar.model$theta==theta & 
                     bivar.model$SF==SF & 
                     bivar.model$rep==rep,]$wmean <- NA
    }
    
}
}
}
}

# Bind area number: 5 groups - E_i [0-20],[20-40],[40-60],[60-80],[80-100]
n <- 5
breaks <- (0:n)/n
bivar.model <- cbind(bivar.model,area=cut(Eprostata,quantile(Eprostata,breaks),include.lowest=T,label=F))

# Save/Load table
# save(bivar.model,file="./bivar_table.RData")
load("bivar_table.RData")

# Plots
Sys.time()
system.time( for( simu.name in simuList) {
for( theta in v.theta) {
for( SF in v.SF) {
  print(paste("> Simu=",simu.name),quote=F)
  print(paste("> Theta=",theta),quote=F)
  print(paste("> SF=",SF),quote=F)
  
  simu <- as.numeric(substr(simu.name,5,5))
  theta.str = gsub("[.]","",as.character(theta))
  scen.res <- bivar.model[bivar.model$simu==simu & 
                bivar.model$theta==theta &
                bivar.model$SF==SF ,]
  
  scen.res.list <- list( bym=scen.res[,c("sel","bym")],
                   crossval=scen.res[,c("sel","crossval")],
                   mean=scen.res[,c("sel","mean")],
                   weigthed_mean=scen.res[,c("sel","wmean")],
                   max=scen.res[,c("sel","max")])
  
  scen.res.list <- lapply(scen.res.list, FUN=function(x) { names(x) <- c("grp","pred"); return(x) })
  
  print(">> ROC curve ...",quote=F)
  title.str <- paste("ROC curves.\nsimu=",simu,",theta=",theta,",SF=",SF)
  roc.plot <- rocplot.multiple(test.data.list=scen.res.list, title=title.str)
  ggsave(plot=roc.plot,file=paste0("./img/roc/roc_",simu.name,"_",theta.str,"_",SF,".png"),width=12,height=10)

  print(">> Scatterplot ...",quote=F)
  title.str <- paste("Scatterplot.\nsimu=",simu,",theta=",theta,",SF=",SF)
  bv.plot <- bivar_scatter_plot(scen.res,title.str,xvar="bym",yvar="crossval")
  ggsave(plot=bv.plot,file=paste0("./img/scat/scat_",simu.name,"_",theta.str,"_",SF,".png"),width=12,height=10)
} # END SF
} # END theta
} # END simu
)
Sys.time()


# Tables
threshold=0.8
fpr.tpr.table <- data.frame()
for( simu in simuList ) {
for( theta in v.theta ) {
for( SF in v.SF ) {
for( rep in 1:10 ) {
  print(paste("> Simu=",simu),quote=F)
  print(paste("> Theta=",theta),quote=F)
  print(paste("> SF=",SF),quote=F)
  
  temp <- bivar.model[ bivar.model$simu==as.numeric( substr(simu,5,5) ) & 
                         bivar.model$theta==theta & 
                         bivar.model$SF==SF,]
  temp <- na.omit(temp)
  
  if( simu=="simu1") {
    # FPR and TPR calculation for each area
    for( area in 1:5 ) {
     temp2 <- temp[temp$area==area,] 
     fpr.tpr.table <- rbind( fpr.tpr.table, 
                             cbind( data.frame( rbind( roc.curve(obs=temp2$sel,pred=temp2$bym,th=threshold),
                                                       roc.curve(obs=temp2$sel,pred=temp2$crossval,th=threshold),
                                                       roc.curve(obs=temp2$sel,pred=temp2$mean,th=threshold),
                                                       roc.curve(obs=temp2$sel,pred=temp2$wmean,th=threshold),
                                                       roc.curve(obs=temp2$sel,pred=temp2$max,th=threshold)), model=c("bym","croosval","mean","wmean","max"), 
                                                simu=as.numeric( substr(simu,5,5) ),rep=rep,theta=theta, SF= SF, area=area ) ) )
     
    }
  }else{
    fpr.tpr.table <- rbind( fpr.tpr.table, 
                            cbind( data.frame( rbind( roc.curve(obs=temp$sel,pred=temp$bym,th=threshold),
                                                      roc.curve(obs=temp$sel,pred=temp$crossval,th=threshold),
                                                      roc.curve(obs=temp$sel,pred=temp$mean,th=threshold),
                                                      roc.curve(obs=temp$sel,pred=temp$wmean,th=threshold),
                                                      roc.curve(obs=temp$sel,pred=temp$max,th=threshold) ), model=c("bym","croosval","mean","wmean","max"), 
                                               simu=as.numeric( substr(simu,5,5) ),rep=rep,theta=theta, SF= SF, area=0 ) ) )
    
  }
  

}
}
}
}

# Save/Load fpr.tpr.table
# save(fpr.tpr.table,file="./fpr_tpr_table.RData")
load("./fpr_tpr_table.RData")

# specificity ( 1-FPR )
dcast(fpr.tpr.table,model+simu+area~theta+SF,mean,value.var="FPR")

# sensitivity
dcast(fpr.tpr.table,model+simu+area~theta+SF,mean,value.var="TPR")

# AUC Table
auc.table <- data.frame()
for( simu in simuList ) {
for( theta in v.theta ) {
for( SF in v.SF ) {
      print(paste("> Simu=",simu),quote=F)
      print(paste("> Theta=",theta),quote=F)
      print(paste("> SF=",SF),quote=F)
      
      for( rep in 1:10 ) {
        print(paste(">> rep=",rep),quote=F)
        
        temp <- bivar.model[ bivar.model$simu==as.numeric( substr(simu,5,5) ) & 
                               bivar.model$theta==theta & 
                               bivar.model$SF==SF,]
        temp <- na.omit(temp)
        
        if( simu=="simu1") {
          # FPR and TPR calculation for each area
          for( area in 1:5 ) {
            temp2 <- temp[temp$area==area,] 
            auc.table <- rbind( auc.table, 
                                    cbind( data.frame( auc=rbind( auc(response=temp2$sel,predictor=temp2$bym),
                                                              auc(response=temp2$sel,predictor=temp2$crossval),
                                                              auc(response=temp2$sel,predictor=temp2$mean),
                                                              auc(response=temp2$sel,predictor=temp2$wmean),
                                                              auc(response=temp2$sel,predictor=temp2$max) ), model=c("bym","crossval","mean","wmean","max"),
                                                       simu=as.numeric( substr(simu,5,5) ),rep=rep,theta=theta, SF= SF, area=area ) ) )
            
          }
        }else{
          auc.table <- rbind( auc.table, 
                                  cbind( data.frame( auc=rbind( auc(response=temp$sel,predictor=temp$bym),
                                                            auc(response=temp$sel,predictor=temp$crossval),
                                                            auc(response=temp$sel,predictor=temp$mean),
                                                            auc(response=temp$sel,predictor=temp$wmean),
                                                            auc(response=temp$sel,predictor=temp$max) ), model=c("bym","crossval","mean","wmean","max"), 
                                                     simu=as.numeric( substr(simu,5,5) ),rep=rep,theta=theta, SF= SF, area=0 ) ) )
          
        }
      }
}
}
}

# Save/Load auc.table
# save(auc.table,file="./auc_table.RData")
load("./auc_table.RData")

dcast(auc.table,model+simu+area~theta+SF,mean,value.var="auc")
