# compare_methods.R
################################
# Author: CNM 
# Date: 18/07/2014
# Project: TFM
# Description: Compare SENSIBILITY and SPECIFICITY values from BYM method and CROSS-VALIDATION method.

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
library(ggplot2)

# Load data
load("iwe_table.RData")
load("spec_sens_table.RData")

compare.table <- rbind(cbind(bym.table,method="BYM"),cbind(iwe.table,method="XVAL") )
compare.table$SF <- paste0( "SF=",compare.table$SF )
compare.table$SF = factor(compare.table$SF, levels=c('SF=1','SF=2','SF=4','SF=10'))
compare.table$simu <- paste( "simu",compare.table$simu )


# Plots
# Spec. = 1-FPR
ggplot(compare.table,aes(y=1-FPR,x=as.factor(theta))) + geom_boxplot(aes(fill=as.factor(method))) + 
  facet_grid(simu~SF) + scale_x_discrete(name=expression(theta)) +
  ggtitle("BYM vs. CROSS VALIDATION method \nSpec. (1-FPR)") +  scale_fill_discrete(name="Method",breaks=c("BYM","XVAL"),labels=c("Besag-York-Mollié","Cross Validation"))


ggplot(compare.table,aes(y=TPR,x=as.factor(theta))) + geom_boxplot(aes(fill=as.factor(method))) + 
  facet_grid(simu~SF) + scale_x_discrete(name=expression(theta)) +
  ggtitle("BYM vs. CROSS VALIDATION method \nTPR (sens.)")  +  scale_fill_discrete(name="Method",breaks=c("BYM","XVAL"),labels=c("Besag-York-Mollié","Cross Validation"))
