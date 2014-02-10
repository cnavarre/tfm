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