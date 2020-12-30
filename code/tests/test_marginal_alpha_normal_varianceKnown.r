source("../pooling_aux.r")

m0 <- -1 ; v0 <- 1^2
m1 <- 2 ; v1 <- 1^2

mv <- c(m0, m1)
vv <- c(v0, v1)
K <- length(mv)

## Observed data
n <- 200
sigmasq <- .02
truemu <- 2
y <- rnorm(n, mean = truemu, sd  = sqrt(sigmasq))
# X <- c(1, 1)/10
X <- c(1.2, 1.2)#/10
M <- 10000
####### Hierarchical priors
alpha.MC.dirichlet <- LearnBayes::rdirichlet(M, X)
alpha.MC.logisticNormal <- rlogisticnorm(N = M,
                                         m = digamma(X)-digamma(X[K]),
                                         Sigma = constructSigma(X))

####### Hierarchical posteriors
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

normaldata.stan <-  list(Y = y, tau = sqrt(sigmasq),
                       X = X, N = n, K = K,
                       mu = mv, sigma_sq = vv)

## Dirichlet
compiled.dirichlet <- stan_model("../stan/posterior_fixedVariance_normal_Dirichlet_pooled.stan")
dirichlet.posterior <- sampling(compiled.dirichlet, data = normaldata.stan, 
                                iter = 5000,
                                control = list(adapt_delta = .99, max_treedepth = 15))
print(dirichlet.posterior, pars = c("theta", "alpha"))

check_hmc_diagnostics(dirichlet.posterior)

mu.dirichlet <- extract(dirichlet.posterior, 'theta')$theta
alphas.dirichlet <- extract(dirichlet.posterior, 'alpha')$alpha
post.alpha.cred.dirichlet <- apply(alphas.dirichlet, 2, quantile, probs = c(.025, .975))
normalPars.dirichlet <- extract(dirichlet.posterior, c('mustar', 'vstar'))
mv.dirichlet <- unlist( lapply(normalPars.dirichlet, mean) ) 
  
## Logistic normal
compiled.logisticNormal <- stan_model("../stan/posterior_fixedVariance_normal_logisticNormal_pooled.stan")

normaldata.stan$means <- digamma(X)-digamma(X[K])
normaldata.stan$Sigma <- constructSigma(X)

logisticNormal.posterior <- sampling(compiled.logisticNormal, data = normaldata.stan,
                                     iter = 5000,
                                     control = list(adapt_delta = .99, max_treedepth = 15))
print(logisticNormal.posterior, pars = c("theta", "alpha"))

check_hmc_diagnostics(logisticNormal.posterior)

mu.logisticNormal <- extract(logisticNormal.posterior, 'theta')$theta
alphas.logisticNormal <- extract(logisticNormal.posterior, 'alpha')$alpha
post.alpha.cred.logisticNormal <- apply(alphas.logisticNormal, 2, quantile, probs = c(.025, .975))
normalPars.logisticNormal <- extract(logisticNormal.posterior, c('mustar', 'vstar'))
mv.logisticNormal <- unlist( lapply(normalPars.logisticNormal, mean) )


#### Testing a simple Monte Carlo approximation of the marginal posterior mean of alpha
logkernel <- function(alpha, data){
  y <- data$y
  sigmasq <- data$s2
  ppars <- pool_par_gauss(alpha, m = mv, v = vv)
  m0 <- ppars[1]
  v0 <- ppars[2]^2
  ###
  xbar <- mean(y)
  s2 <-  sum((y-xbar)^2)
  vtilde <- n*v0 + sigmasq 
  l1 <- -s2/(2*sigmasq)
  l2 <- - (n*(xbar-m0)^2)/(2*vtilde)
  l3 <- -0.5 * ( log(vtilde) + (n-1) * log(2*pi*sigmasq) + log(2*pi))
  ans <- l1 + l2 + l3
  return(ans)
}

lZ.dirichlet <- matrixStats::logSumExp(apply(alpha.MC.dirichlet, 1, logkernel, data = list(y = y, s2 = sigmasq)))-log(M) 
lZ.logisticNormal <- matrixStats::logSumExp(apply(alpha.MC.logisticNormal, 1, logkernel, data = list(y = y, s2 = sigmasq)))-log(M) 

lZ.dirichlet
bridgesampling::bridge_sampler(dirichlet.posterior)

lZ.logisticNormal
bridgesampling::bridge_sampler(logisticNormal.posterior)

phi <- function(alpha, data, lZ){
  y <- data$y
  sigmasq <- data$s2
  ppars <- pool_par_gauss(alpha, m = mv, v = vv)
  m0 <- ppars[1]
  v0 <- ppars[2]^2
  ###
  xbar <- mean(y)
  s2 <-  sum((y-xbar)^2)
  vtilde <- n*v0 + sigmasq 
  l1 <- -s2/(2*sigmasq)
  l2 <- - (n*(xbar-m0)^2)/(2*vtilde)
  l3 <- -0.5 * ( log(vtilde) + (n-1) * log(2*pi*sigmasq) + log(2*pi))
  ans <- log(alpha) + l1 + l2 + l3 - lZ
  ans <- exp(ans)
  return(ans)
}

(post.phi.dirichlet <- rowMeans(apply(alpha.MC.dirichlet, 1, phi, data = list(y = y, s2 = sigmasq), lZ  = lZ.dirichlet)) )
(post.phi.logisticNormal <- rowMeans(apply(alpha.MC.logisticNormal, 1, phi, data = list(y = y, s2 = sigmasq), lZ = lZ.logisticNormal)) )

colMeans(alphas.dirichlet)
colMeans(alphas.logisticNormal)

###########
# a closer look at the marginal posterior for alpha
# Prior
a0 <- X[1]
b0 <- X[2]
hist(alpha.MC.dirichlet[, 1], probability = TRUE)
curve(dbeta(x, a0, b0), lwd = 2, add = TRUE)

## Posterior
marginal_alpha <- function(x, log = FALSE){
  ppars <- pool_par_gauss(c(x, 1-x), m = mv, v = vv)
  ans <- normal_mean_marg_like(y = y, sigma = sqrt(sigmasq), m = ppars[1], v = ppars[2]^2, log = TRUE)  + dbeta(x, a0, b0, log = TRUE)
  if(!log) ans <- exp(ans)
  return(ans)
}
marginal_alpha <- Vectorize(marginal_alpha)

Z <- integrate(marginal_alpha, 0, 1)$value

marginal_alpha_norm <- function(x, log = FALSE){
  ans <- marginal_alpha(x, log = TRUE) - log(Z)
  if(!log) ans <- exp(ans)
  return(ans)
}
marginal_alpha_norm <- Vectorize(marginal_alpha_norm)

ZZ <- integrate(marginal_alpha_norm, 0, 1)$value
marginal_alpha_norm2 <- function(x, log = FALSE){
  ans <- marginal_alpha_norm(x, log = TRUE) - log(ZZ)
  if(!log) ans <- exp(ans)
  return(ans)
}
marginal_alpha_norm2 <- Vectorize(marginal_alpha_norm2)
integrate(function(x) marginal_alpha_norm2(x), 0 , 1)

alpha.posterior.samples <- alphas.dirichlet
hist(alpha.posterior.samples, probability = TRUE, main = paste("n=", n, " true mean:", truemu, " m0:", m0))
curve(marginal_alpha_norm2, lwd = 2, add = TRUE)
curve(dbeta(x, a0, b0), lwd = 2, lty = 2,  add = TRUE)

####
# E[alpha_0]
colMeans(alphas.dirichlet)[1]
post.phi.dirichlet[1]
(mu.quad <- integrate(function(x) x* marginal_alpha_norm2(x), 0 , 1) )
# SD(alpha_0)
ex2 <- integrate(function(x) x*x* marginal_alpha_norm2(x), 0 , 1)$value
sqrt(ex2 - mu.quad$value^2)
sd(alphas.dirichlet[, 1])