source("../pooling_aux.r")
# #########
m0 <- 0 ; v0 <- .25
m1 <- 0; v1 <- .75
m2 <- 0; v2 <- 1
m3 <- 0; v3 <- 1
m4 <- 0; v4 <- .5
mus <- c(m0, m1, m2, m3, m4)
sigmas <-  c(v0, v1, v2, v3, v4) 
K <- length(sigmas)
set.seed(666)
smp <- sample(1:100, K, replace = TRUE) 
( alphas <- smp/sum(smp) )
ws <- alphas/sigmas
mstar <- sum(ws*mus)/sum(ws)
vstar <-  1/sum(ws)
########
# Preliminaries
gY <- function(y, mu, v){
  omega <- c(-sqrt(y), sqrt(y))
  dens <- sapply(omega, function(x) dnorm(x, mu, sqrt(v))/abs(2 * x))
  return(sum(dens))
}
gY <- Vectorize(gY)

## (a) Pool-then-induce
piY <- function(y, vst){
  gY(y, mstar, vst)
}
piY <- Vectorize(piY)
## (b) Induce-then-pool
Ds <- list(
  f0 = function(x){gY(x, m0, v = v0)},
  f1 = function(x){gY(x, m1, v = v1)},
  f2 = function(x){gY(x, m2, v = v2)},
  f3 = function(x){gY(x, m3, v = v3)},
  f4 = function(x){gY(x, m4, v = v4)}
) # list with the densities

xlwr <- 0
xupr <- 5
curve(piY(x, vstar),
      cex.lab = 1.5, cex.axis = 1.5,
      xlwr, xupr, main = "", xlab = expression(Y),
      ylab = "Density", lwd = 3)
curve(dpoolnorm.positive(x, D = Ds, alpha = alphas), xlwr, xupr, lwd = 3,
      lty = 2, col = "grey50", add = TRUE)
legend(x = "topright", legend = c(expression(pi), expression(pi^"|")), 
       col = c("black", "grey50"), lty = 1:2, bty = "n", lwd = 3, cex = 2)

t(
  sapply(seq(0, 10, length.out = 10), function(x){
    return(
      c(
        dpoolnorm.positive(x, D = Ds, alpha = alphas),
        piY(x, vst = vstar)
      )
    )
  })
)
