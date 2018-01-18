source("pooling_aux.r")
# #########
# m0 <- 6 ; v0 <- .25
# m1 <- 5.5; v1 <- .25
# m2 <- 6; v2 <- 1
# m3 <- 10; v3 <- 1
# m4 <- 5; v4 <- .5
m0 <- 6E3; v0 <- 5E2^2
m1 <- 5.5E3; v1 <- 2E3^2
m2 <- 6E3; v2 <- 1E3^2
m3 <- 1E4; v3 <- 1E3^2
m4 <- 5E3; v4 <- 5E3^2
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
Ds <- list(
  f0 = function(x){dnorm(x, mean = m0, sd = sqrt(v0))},
  f1 = function(x){dnorm(x, mean = m1, sd = sqrt(v1))},
  f2 = function(x){dnorm(x, mean = m2, sd = sqrt(v2))},
  f3 = function(x){dnorm(x, mean = m3, sd = sqrt(v3))},
  f4 = function(x){dnorm(x, mean = m4, sd = sqrt(v4))}
) # list with the densities

xlwr <- 0
xupr <- 10E3 
curve(dnorm(x, mean = mstar, sd = sqrt(vstar)),
      xlwr, xupr, main = "Pooled distributions", xlab = expression(theta),
      ylab = "Density", lwd = 2)
curve(dpoolnorm(x, D = Ds, alpha = alphas), xlwr, xupr, lwd = 2,
      lty = 2, add = TRUE, col = 2)
legend(x = "topright", bty = "n", legend = c("Analytic", "Computational"),
       col = 1:2, lty = 1:2, lwd = 2) ## Maths works!