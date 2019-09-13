alpha0 <- .7
alpha1 <- 1-alpha0
alphas <- c(alpha0, alpha1)

mU0 <-  -.80 ; vU0 <- .40
mU1 <- .5; vU1 <- .05

mV0 <- -1.6; vV0 <- 0.024
mV1 <- -1.25; vV1 <- 0.4

M <- 1e6
U0 <- rlnorm(M, m = mU0, sd = sqrt(vU0))
U1 <- rlnorm(M, m = mU1, sd = sqrt(vU1))

V0 <- rlnorm(M, m = mV0, sd = sqrt(vV0))
V1 <- rlnorm(M, m = mV1, sd = sqrt(vV1))

Z0 <- U0/V0

Z1 <- U1/V1

hist(Z0)
hist(Z1)

muU <- c(mU0, mU1)
muV <- c(mV0, mV1)
sigmaU <- c(vU0, vU1)
sigmaV <- c(vV0, vV1)

###########
## Pool-then-induce

wsU <- alphas/sigmaU
muUstar <- sum(wsU*muU)/sum(wsU)
vUstar <-  1/sum(wsU)
#
wsV <- alphas/sigmaV
muVstar <- sum(wsV*muV)/sum(wsV)
vVstar <-  1/sum(wsV)

muZstar <- muUstar - muVstar  
vZstar <- vUstar + vVstar

## Induce then pool
wsZ <- alphas/(sigmaU + sigmaV)
muZstarstar <- sum(wsZ * muU)/sum(wsZ) - sum(wsZ * muV)/sum(wsZ)
vZstarstar <- 1/sum(wsZ)

curve(dlnorm(x, muZstar, sqrt(vZstar)) , 0, 15, ylim = c(0, .25),
      ylab = "Density", xlab = expression(Z), cex.lab = 1.5, cex.axis = 1.5)
curve(dlnorm(x, muZstarstar, sqrt(vZstarstar)), lty = 2, add = TRUE)
legend(x = "topright", legend = c(expression(pi), expression(pi^"|")), lty = 1:2, bty = "n", cex = 2)

