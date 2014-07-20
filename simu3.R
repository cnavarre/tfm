# DEFINITION: 
# Situation with high heterogeneity compising NOC cluster of 1% of expected values.
# Observations are simulated from multinomial of size N=[sum(E_i)] and p_i {prop} E_i*theta_i
# INPUT PARAMETERS:
#  v.exp <- vector of expected values for background observed values
#  all.nb <- nb object
#  NOC -> number of clusters
#  theta -> value of theta, modifier of expeted value for clusters
#  SF -> scale factor for expected values for all region.
#  pct -> minimun percentage of expected values in cluster.
# OUTPUT VALUES:
#  exp <- vector of expected values used in Multinomial
simu3 <- function( v.exp,all.nb,theta,NOC=20,SF=1,pct=0.01 ) {
  require(spdep);require("maptools")
  
  if( pct<0 | pct>1) stop('Not valid percentage value for cluster')
  
  # Factor de escala para los valores esperados
  v.exp <- v.exp * SF
  v.exp.base <- v.exp
  
  tot.exp <- sum(v.exp) # Total amount of expected values
  N <- floor(tot.exp)
  
  # IDs
  ID <- names(v.exp)

  # Vector de areas incluidas en los clusters
  sel.idx <- list()
  
  while( length(sel.idx) < NOC ){
    
    # New cluster
    clus.idx <- sample( ID[card(all.nb)>0], size=1 )
    
    cont <- TRUE  
    while( cont ) {
      list.nb <- cluster.nb( which( (ID %in% clus.idx) ), all.nb ) # Neighbours of last added
      
      # Spatial units to add into cluster
      obj2add.nb <- list.nb[cumsum(( v.exp[list.nb] + sum(v.exp[clus.idx]) ) / tot.exp) < pct]
      clus.idx <- c(clus.idx,ID[obj2add.nb])
      
      # Check if cluster reach pct% of total amount of expected values
      cont <- ( length(obj2add.nb) == length(list.nb) )
    }
  
    # Check if new cluster (or his frontier) is nonoverlapping with previous clusters
    # in this case: add to selected areas to increase the risk
    clus.front.idx <- c(clus.idx, ID[list.nb] ) # cluster and frontier's cluster
    if( !any( clus.front.idx %in% unlist(sel.idx) ) ) 
      sel.idx[[length(sel.idx)+1]] <- clus.idx      
  }
  
  v.obs <- vector()
  
  # Asignamos el valor multiplicado por theta para cada punto seleccionado
  v.exp[ ID %in% unlist(sel.idx) ] <- v.exp[ ID %in% unlist(sel.idx) ] * theta
# v.exp <- v.exp / sum(v.exp) # normalization not necessary for generation
  
  # Creamos el vector de casos observados
  v.obs <- as.vector(rmultinom(1,size=N,prob=v.exp))
  names(v.obs) <- ID
  
  # Devolvemos una lista con el valor esperado,valor observado y el identificador
  sel.idx <- lapply(sel.idx,FUN=function(x){ as.numeric(paste("46250",x,sep="")) })
  return(list(exp=v.exp,obs=v.obs,sel=sel.idx,theta=theta,SF=SF,exp.base=v.exp.base))
}


