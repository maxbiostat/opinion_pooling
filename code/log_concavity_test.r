#######
## A certain result in the paper says that if the f_i are log-concave, pi() will also be.
## This is a successful test of that fact.
#######
source("pooling_aux.r")
###
Ds <- list(
  f0 = function(x){dnorm(x, mean = 0, sd = 5)},
  f1 = function(x){dnorm(x, 10, 2)},
  f2 = function(x){dexp(x, rate = 1/5)},
  f3 = function(x){dlogis(x, location = -15, scale = 2)}
) # list with the densities
collection.dens <- function(x, alphas){
  return(dpoolnorm(x, D = Ds, alphas = alphas))
}
lapply(Ds, function(f) curve(f(x), -20, 20))
Alphas <- rep(1/4, 4)
poolFun <- function(y) collection.dens(y, alphas = Alphas)
curve(poolFun, -20, 20, lwd = 2, col = "green")