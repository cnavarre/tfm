#' @description Creates simulation scenario type.2  
#' \code{simu2} creates an scenario type.2
#' @details
#' 
#' @param v.exp a vector of expected values for background observed values
#' @param all.nb a nb object
#' @param theta modifier of expeted value for cluster elements
#' @param SF scale factor for expected values for all region.
#' @param pct is minimun percentage of expected values in cluster.
#' @return exp modified vector with expected values
#' @return obs vector of observations
#' @return sel vector of indexes of areas selected into the cluster
#' @example 
#' simu2(v.exp=Eprostata,all.nb=vlc.nb,theta=2,SF=1,pct=0.01)

simu2 <- function( v.exp,all.nb,theta,SF=1,pct=0.01 ) {
  require(spdep);require("maptools")
  
  # if( pct<0 | pct>1) stop("Not valid pct. value for cluster")
  
  # Factor de escala para los valores esperados
  v.exp <- v.exp * SF
  tot.exp <- sum(v.exp) # Total amount of expected values
  
  # IDs
  ID <- names(v.exp)

  # Seleccionamos un barrio al azar con conexiones
  sel.idx <- vector()
  sel.idx <- sample( ID[card(all.nb)>0], size=1 )
  
  cont <- TRUE  
  while( cont ) {
    list.nb <- cluster.nb( which( (ID %in% sel.idx) ), all.nb ) # Neighbours of last added
    
    # To be added to cluster
    obj2add.nb <- list.nb[cumsum(( v.exp[list.nb] + sum(v.exp[sel.idx]) ) / tot.exp) < pct]
    sel.idx <- c(sel.idx,ID[obj2add.nb])
    
    cont <- ( length(obj2add.nb) == length(list.nb) )
  }

  v.obs <- vector()
  
  # Asignamos el valor multiplicado por theta para cada punto seleccionado
  v.exp[ ID %in% sel.idx ] <- v.exp[ ID %in% sel.idx ] * theta
  
  # Creamos el vector de observaciones
  v.obs <- sapply(v.exp,FUN=function(lambda){rpois(1,lambda)})
  names(v.obs) <- ID
  
  # Devolvemos una lista con el valor esperado,valor observado y el identificador
  sel.idx <- as.numeric(paste("46250",sel.idx,sep=""))
  return(list(exp=v.exp,obs=v.obs,sel=sel.idx))
}
