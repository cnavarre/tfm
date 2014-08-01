# Ranking vector values by percentil
perc.rank <- function(x) trunc(rank(x))/length(x)

# Function that returns a cluster neighbours list
cluster.nb <- function( cluster,all.nb ){ 
  # Get all neighbours
  list.nb <- unique(unlist(lapply(cluster,FUN=function(idx){ all.nb[[idx]] } )))
  # Remove 0 values -> members of cluster without nb's.
  list.nb <- list.nb[list.nb!=0] 
  # Remove nb's in cluster
  list.nb <- list.nb[!(list.nb %in% cluster)]
  
  return(list.nb)
}

destring <- function(x) {
  ## convert factor to strings
  if (is.character(x)) {
    as.numeric(x)
  } else if (is.factor(x)) {
    as.numeric(levels(x))[x]
  } else if (is.numeric(x)) {
    x
  } else {
    stop("could not convert to numeric")
  }
}


.eval <- function(evaltext,envir=sys.frame()) {
  ## evaluate a string as R code
  eval(parse(text=evaltext), envir=envir)
}

## trim white space/tabs
## this is marek's version
trim<-function(s) gsub("^[[:space:]]+|[[:space:]]+$","",s)

## auto-install packages 
## if not found
libra <- function(x) { 
  if (!base::require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE) ; 
    base::require(x, character.only = TRUE)
  } 
}

## remove if exist
rem <- function( varList ) {
  for( varName in varList )   
    if( exists(varName) ) rm(list=varName,envir=parent.frame())
}

## Classifing areas through SMR posterior matrix
## High risk (1): if P( \theta_i > R0 ) > 1
## Low Risk (0): otherwise
decision.rule <- function( smr.m , k=0.8, R0=1 ) {
  prob <- apply(smr.m > R0,2,sum) / apply(smr.m,2,length)
  risk.level <- factor(as.numeric(prob > k),levels=c("0","1"))
  return(risk=risk.level)
}


## Create a NULL matrix of dimension nrow x ncol
create.matrix <- function(nrow, ncol) {
  x <- matrix()
  length(x) <- nrow * ncol
  dim(x) <- c(nrow,ncol)
  x
}


# Colnames generator
col_name <- function( name, number ) {
  paste0( name, "[", number, "]" )
}

# ROC curve table
roc.curve=function(obs,pred,th,print=FALSE){
  Y=obs
  Ps=(pred>th)*1
  FP=sum((Ps==1)*(Y==0))/sum(Y==0)
  TP=sum((Ps==1)*(Y==1))/sum(Y==1)
  if(print){
   print(table(Observed=Y,Predicted=Ps))
  }
  vect=c(FP,TP)
  names(vect)=c("FPR","TPR")
  return(vect)
}


# Plot scenario
plot.scenario <- function( scenFile="default.scen",data=NULL, map=Carto,v.exp=Eprostata, plot=TRUE, save=FALSE ) {
  require(ggplot2)
  require(raster)
  require(rgeos)
  
  # Read scenario data
  if(is.null(data) ){
    data <- dget( scenFile )
  }
  
  sel.idx <- unlist( lapply( unlist(data$sel), FUN=function(x) which(names(data$obs)==substr(x,6,10)) ) )
  
  xlim = c(720678,735386)
  ylim = c(4362066,4383252)
  
  # Create the clipping polygon
  CP <- as(extent(xlim[1], xlim[2], ylim[1], ylim[2]), "SpatialPolygons")
  proj4string(CP) <- CRS(proj4string(Carto))
  
  # Clip the map
  map.clip <- gIntersection(map, CP, byid=TRUE,id=rownames(map@data))
  
  # Create the polygons
  map.ggplot <- fortify(map.clip)
  # Merge the Cases value to the polygons
  # map.ggplot <- cbind(map.ggplot, exp = map[map.ggplot$id, "exp"])
  map.ggplot$sel <- map.ggplot$id %in% unlist(data$sel)
  map.plot <- ggplot(data = map.ggplot) + geom_polygon(aes(x = long, y = lat, fill = sel, group = group)) + 
            coord_equal() + scale_fill_brewer("Selected area", type="qual",palette=2) + 
            ggtitle(paste(basename(scenFile),"\n",round(sum(v.exp[sel.idx])/sum(v.exp)*100,2),"% expected cases over total cases in cluster.")) + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
  if(plot) map.plot
  if(save) { suppressMessages( ggsave(filename=paste0(getwd(),"/img/",sapply(strsplit(basename(scenFile),"\\."), function(x) paste(x[1:(length(x)-1)], collapse=".")),".png"), map.plot )  ) }
  return(map.plot)
}



# Bivariate scatter plot with histograms
bivar_scatter_plot2 <- function( data ) {
  require(gridExtra)
  require(ggplot2)
  
  hist_top <- ggplot(data) + geom_histogram( aes(x=log(x)) )
  empty <- ggplot() + geom_point(aes(1,1), colour="white") + theme_bw() +
    theme(axis.ticks=element_blank(), 
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),           
          axis.title.x=element_blank(), axis.title.y=element_blank())
  
  hist_right <- ggplot(data) + geom_histogram(aes(x=-log(y))) + coord_flip()
  
  scatter <- ggplot(data) + geom_point( aes( x=log(x), y=-log(y), colour=as.factor(sel) ) ) + scale_colour_discrete(guide = FALSE)
  
  grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}

# Another interesting plot for bivariate distribution
bivar_scatter_plot <- function( data, title,xvar=x,yvar=y ) {
  require(ggplot2)
  
  scatter <- ggplot(data,aes_string(x=xvar,y=yvar)) + geom_point(aes(color=as.factor(sel)),alpha=0.7 ) + scale_color_manual( values=c("orange","purple")) +
    scale_x_continuous() + 
    scale_y_continuous() + 
    geom_rug(col=rgb(255/255,165/255,0,alpha=.2)) + ggtitle(title) + theme(legend.position="none")
  return(scatter)
}



# ROC data
rocdata <- function(grp, pred){
  # Produces x and y co-ordinates for ROC curve plot
  # Arguments: grp - labels classifying subject status
  #            pred - values of each observation
  # Output: List with 2 components:
  #         roc = data.frame with x and y co-ordinates of plot
  #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval
  
  grp <- as.factor(grp)
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  } 
  
  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }
  
  cut <- unique(pred)
  tp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
  fn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
  fp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
  tn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc = data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]
  
  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
  
  pos <- pred[grp == levels(grp)[2]]
  neg <- pred[grp == levels(grp)[1]]
  q1 <- auc/(2-auc)
  q2 <- (2*auc^2)/(1+auc)
  se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) -1)*(q1 - auc^2)) + ((length(neg) -1)*(q2 - auc^2)))/(length(pos)*length(neg)))
  ci.upper <- auc + (se.auc * 0.96)
  ci.lower <- auc - (se.auc * 0.96)
  
  se.auc.null <- sqrt((1 + length(pos) + length(neg))/(12*length(pos)*length(neg)))
  z <- (auc - 0.5)/se.auc.null
  p <- 2*pnorm(-abs(z))
  
  stats <- data.frame (auc = auc,
                       p.value = p,
                       ci.upper = ci.upper,
                       ci.lower = ci.lower
  )
  
  return (list(roc = roc, stats = stats))
}



# Function to plot a single ROC curve
rocplot.single <- function(grp, pred, title = "ROC Plot", p.value = FALSE){
  require(ggplot2)
  plotdata <- rocdata(grp, pred)
  
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
  }
  
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
    geom_line(aes(colour = "")) +
    geom_abline (intercept = 0, slope = 1) +
    theme_bw() +
    scale_x_continuous("False Positive Rate (1-Specificity)") +
    scale_y_continuous("True Positive Rate (Sensitivity)") +
    scale_colour_manual(labels = annotation, values = "#000000") +
    ggtitle(title) +
    theme( plot.title = element_text(face="bold", size=14), 
         axis.title.x = element_text(face="bold", size=12),
         axis.title.y = element_text(face="bold", size=12, angle=90),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.justification=c(1,0), 
         legend.position=c(1,0),
         legend.title=element_blank(),
         legend.key = element_blank()
    )
  return(p)
}


# Plots multiple ROC curves
rocplot.multiple <- function(test.data.list, groupName = "grp", predName = "pred", title = "ROC Plot", p.value = TRUE){
  require(plyr)
  require(ggplot2)
  plotdata <- llply(test.data.list, function(x) with(x, rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
  plotdata <- list(roc = ldply(plotdata, function(x) x$roc),
                   stats = ldply(plotdata, function(x) x$stats)
  )
  
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
  }
  
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
    geom_line(aes(colour = .id)) +
    geom_abline (intercept = 0, slope = 1) +
    theme_bw() +
    scale_x_continuous("False Positive Rate (1-Specificity)") +
    scale_y_continuous("True Positive Rate (Sensitivity)") +
    scale_colour_brewer(palette="Set1", breaks = names(test.data.list), labels = paste(names(test.data.list), ": ", annotation, sep = "")) +
    ggtitle(title) +
    theme(plot.title = element_text(face="bold", size=14), 
         axis.title.x = element_text(face="bold", size=12),
         axis.title.y = element_text(face="bold", size=12, angle=90),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.justification=c(1,0), 
         legend.position=c(1,0),
         legend.title=element_blank(),
         legend.key = element_blank()
    )
  return(p)
}

