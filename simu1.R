# Load data
setwd("~/docs/cursos/master_uv_bioestadistica/project/scripts/")
load("./data/input_data.Rdata")

# Load shared functions
setwd("~/GitHub/tfm")
source("./utils.tfm.R")

# DEFINITION: simu1
# Function to generate simulated scenario with n spatial units with E_i*theta as new
# expected values. Final scenario has 'n' single isolated areas with elevated risks.
# Observed number of cases (Y_i) comes from Poisson( E_i*theta_i ), where theta_i = 1 if 
# the area is not in selected indexes vector (sel.idx), otherwise theta_i=theta parameter.
#
# INPUT: 
#  v.exp -> E_i, background expected observations for each geo-unit (i).
#  n -> number of units to increase the expectation value.
#  theta -> proportional constant.
#  SF -> scale factor
#
# OUTPUT:
#  exp -> expected values for each spatial unit
#  obs -> number of observed cases per spatial unit.
#  sel.idx -> selected indexes with new expected value: E^*_i = E_i*theta
simu1 <- function( v.exp,n,theta,SF=1 ) {
  # Factor de escala para los valores esperados
  v.exp <- v.exp * SF
  
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


simu1(Eprostata,5,1.5)

simu1(Eprostata,5,2,SF=2)

simu1(Eprostata,5,3,SF=10)
