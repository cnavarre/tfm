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