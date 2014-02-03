load("./data/input_data.Rdata")

# Definicion del vector de thetas
v.theta = c(1.5,2,3)

# Factor de escala para la los valores esperados
v.SF = 2:10

simu1 <- function( v.exp,n,theta,SF=1 ) {
  # Factor de escala para los valores esperados
  v.exp <- v.exp * SF
  
  # Calculamos los ranking por percentil de cada valor del vector:
  perc.rank <- function(x) trunc(rank(x))/length(x)
  
  df.prost <- data.frame( exp=v.exp )
  df.prost <- within(df.prost, prank <- perc.rank(v.exp))
  df.prost$id <- names(v.exp)

  # Separamos la muestra en 'n' trozos y seleccionamos una muestra de cada trozo
  breaks <- (0:n)/n
  sel.idx <- vector()
  for( i in 2:length(breaks) ) {
    sel.idx <- c(sel.idx,
                 sample( df.prost[df.prost$prank <= breaks[i] & 
                    df.prost$prank > breaks[i-1], "id" ], size=1 ) )
  }

  v.obs <- vector()
  
  # Asignamos el valor multiplicado por theta para cada punto seleccionado
  v.exp[ df.prost$id %in% sel.idx ] <- 
    v.exp[ df.prost$id %in% sel.idx ] * theta
  
  # Creamos el vector de observados mediante una Poisson
  v.obs <- sapply(v.exp,FUN=function(lambda){rpois(1,lambda)})
  names(v.obs) <- df.prost$id
  
  # Devolvemos una lista con el valor esperado (modificado con theta),valor observado y el identificador
  return(list(exp=v.exp,obs=v.obs,sel=sel.idx))
}

simu1(Eprostata,5,v.theta[1])

simu1(Eprostata,5,v.theta[2],SF=2)
