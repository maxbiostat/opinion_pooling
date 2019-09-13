source("pooling_aux.r")
devtools::source_url("https://raw.githubusercontent.com/maxbiostat/CODE/master/R/DISTRIBUTIONS/lognormal_parameter_transformation.r")
###########
pool_par_lognormal <- function(alpha, m, v){
  ## same as pool_par(), but for the Gaussian distribution
  ## Takes in MEANS and VARIANCES and outputs MEAN and SD (for plotting reasons)
  ws <- alpha/v
  vstar <-  1/sum(ws)
  mstar <- sum(ws*m) * vstar
  c(mstar, vstar)
}
compute_pars <- function(alphas, mU, mV, sU, sV){
  ## Pool-then-induce
  PtI.pars.U <- pool_par_lognormal(alphas, mU, sU)
  PtI.pars.V <- pool_par_lognormal(alphas, mV, sV)
  muZstar <- PtI.pars.U[1] - PtI.pars.V[1]  
  vZstar <- PtI.pars.U[2] + PtI.pars.V[2]
  ## Induce then pool
  wsZ <- alphas/(sU + sV)
  muZstarstar <- sum(wsZ * mU)/sum(wsZ) - sum(wsZ * mV)/sum(wsZ)
  vZstarstar <- 1/sum(wsZ)
  return(list(
    muZstar = muZstar,
    vZstar = vZstar,
    muZstarstar = muZstarstar,
    vZstarstar = vZstarstar
  ))
}
KL_lognormals <- function(m1, v1, m2, v2){
  ## m_i is logmean and s_i is logVARIANCE
  1/(2*v2) * {(m1-m2)^2 + v1-v2} + log(v1)-log(v2) 
}
get_transformed_KL <- function(alphas, mU, sU, mV, sV, direction = c("ab", "ba")){
  dist.pars <- compute_pars(alphas, mU = muU, sU = sigmaU, mV = muV, sV = sigmaV)
  ## Directions: (a) = Pool-then-induce (star) and (b) Induce-then-pool (starstar)
  KL <- switch(
    direction,
    "ab" = KL_lognormals(m1 = dist.pars$muZstar, v1 = dist.pars$vZstar, m2 = dist.pars$muZstarstar, v2 = dist.pars$vZstarstar),
    "ba" = KL_lognormals(m1 = dist.pars$muZstarstar, v1 = dist.pars$vZstarstar, m2 = dist.pars$muZstar, v2 = dist.pars$vZstar)
  )
  return(KL)
}
###################################

mU0 <- LogMean(1.5, 0.1)
vU0 <- LogVar(1.5, 0.1)
  
mU1 <- LogMean(1.5, 0.5)
vU1 <- LogVar(1.5, 0.5)
  
mU2 <- LogMean(realMean = 2.0, realSD = 0.5)
vU2 <- LogVar(realMean = 2.0, realSD = 0.5)

mV0 <- LogMean(realMean = 0.5, realSD = 0.1)
vV0 <- LogVar(realMean = 0.5, realSD = 0.1)

mV1 <- LogMean(realMean = 0.5, realSD = 0.2)
vV1 <- LogVar(realMean = 0.5, realSD = 0.2)

mV2 <- LogMean(realMean = 0.7, realSD = 0.01)
vV2 <- LogVar(realMean = 0.7, realSD = 0.01)


muU <- c(mU0, mU1, mU2)
muV <- c(mV0, mV1, mV2)
sigmaU <- c(vU0, vU1, vU2)
sigmaV <- c(vV0, vV1, vV2)

K <- length(sigmaU)

####### Minimise KL in transformed space
optkl_lognorm_a <- function(alphaR){
  get_transformed_KL(alphas = alpha_01(alphaR), mU = muU, sU = sigmaU, mV = muV, sV = sigmaV, direction = "ab")
}
optkl_lognorm_b <- function(alphaR){
  get_transformed_KL(alphas = alpha_01(alphaR), mU = muU, sU = sigmaU, mV = muV, sV = sigmaV, direction = "ba")
}
###############

N <- 1000 ## could increase to, say, 10000 in order to make sure, but it's fine
kl.many.startingPoints <- matrix(rnorm(n = (K-1)*N, mean = 0, sd = 100), ncol = K-1, nrow = N)
many.kls <- lapply(1:N, function(i) {
  optim(kl.many.startingPoints[i, ], optkl_lognorm_a)
})
optimised.kls <- unlist(lapply(many.kls, function(x) x$value))

hist(optimised.kls)
abline(v = optimised.kls[which.min(optimised.kls)], lty = 2, lwd = 2)

( alphaKL.opt <- alpha_01(many.kls[[which.min(optimised.kls)]]$par) )


lognormal.ratio.pars <- compute_pars(alphaKL.opt, mU = muU, sU = sigmaU, mV = muV, sV = sigmaV)

curve(dlnorm(x, lognormal.ratio.pars$muZstar, sqrt(lognormal.ratio.pars$vZstar)) , 0, 10, ylim = c(0, .75),
      ylab = "Density", xlab = expression(Z), cex.lab = 1.5, cex.axis = 1.5)
curve(dlnorm(x, lognormal.ratio.pars$muZstarstar, sqrt(lognormal.ratio.pars$vZstarstar)), lty = 2, add = TRUE)
legend(x = "topright", legend = c(expression(pi), expression(pi^"|")), lty = 1:2, bty = "n", cex = 2)

