load("./data/input_data.Rdata")

# Definicion del vector de thetas
v.theta = c(1.5,2,3)

# Factor de escala para la los valores esperados
v.SF = 2:10

# nb object from spatial polygons.
require(spdep);require("maptools")
vlc.nb <- poly2nb(Carto,snap=0.01)


# DEFINITION: 
# INPUT PARAMETERS:
#  v.exp <- vector of expected values for background observed values
#  all.nb <- nb object
#  NOC -> number of clusters
#  theta -> value of theta, modifier of expeted value for clusters
#  SF -> scale factor for expected values for all region.
#  pct -> minimun percentage of expected values in cluster.
# OUTPUT VALUES:
#  v.exp <- 
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
    n <- length(sel.idx)
    last.added <- sel.idx[n]
    list.nb <- all.nb[[which( (ID %in% last.added) )]] # Neighbours of las added
    list.nb <- list.nb[!(ID[ list.nb ] %in% sel.idx )] # Filter out already added in previous steps
    # To be added to cluster
    obj2add.nb <- list.nb[cumsum(( v.exp[list.nb] + sum(v.exp[sel.idx]) ) / tot.exp) < pct]
    sel.idx <- c(sel.idx,ID[obj2add.nb])
    
    cont <- ( length(obj2add.nb) == length(list.nb) )
  }

  v.obs <- vector()
  
  # Asignamos el valor multiplicado por theta para cada punto seleccionado
  v.exp[ ID %in% sel.idx ] <- v.exp[ ID %in% sel.idx ] * theta
  
  # Creamos el vector de 
  v.obs <- sapply(v.exp,FUN=function(lambda){rpois(1,lambda)})
  names(v.obs) <- ID
  
  # Devolvemos una lista con el valor esperado,valor observado y el identificador
  return(list(exp=v.exp,obs=v.obs,sel=sel.idx))
}


# Testing
simu2(v.exp=Eprostata,all.nb=vlc.nb,theta=v.theta[1],SF=1,pct=0.01)
