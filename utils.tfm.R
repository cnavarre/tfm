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
