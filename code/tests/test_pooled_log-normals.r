source("../pooling_aux.r")
# #########
m0 <- -6 ; v0 <- .25
m1 <- -5.5; v1 <- .25
m2 <- -6; v2 <- 1
m3 <- -10; v3 <- 1
m4 <- -5; v4 <- .5
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
  f0 = function(x){dlnorm(x, meanlog = m0, sdlog = sqrt(v0))},
  f1 = function(x){dlnorm(x, meanlog = m1, sdlog = sqrt(v1))},
  f2 = function(x){dlnorm(x, meanlog = m2, sdlog = sqrt(v2))},
  f3 = function(x){dlnorm(x, meanlog = m3, sdlog = sqrt(v3))},
  f4 = function(x){dlnorm(x, meanlog = m4, sdlog = sqrt(v4))}
) # list with the densities

xlwr <- 0
xupr <- 1E-2 
curve(dlnorm(x, meanlog = mstar, sdlog = sqrt(vstar)),
      xlwr, xupr, main = "Pooled distributions", xlab = expression(theta),
      ylab = "Density", lwd = 2)
curve(dpoolnorm(x, D = Ds, alpha = alphas), xlwr, xupr, lwd = 2,
      lty = 2, add = TRUE, col = 2)
legend(x = "topright", bty = "n", legend = c("Analytic", "Computational"),
       col = 1:2, lty = 1:2, lwd = 2) ## Maths works!