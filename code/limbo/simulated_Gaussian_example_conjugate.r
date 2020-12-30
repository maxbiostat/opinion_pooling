source("pooling_aux.r")
library(rstan)
rstan_options(auto_write = TRUE)

compiled.dirichlet <- stan_model("stan/posterior_conjugate_normal_Dirichlet_pooled.stan")
compiled.logisticNormal <-   stan_model("stan/posterior_conjugate_normal_logisticNormal_pooled.stan")

###################
parameters <- list(
  m = c(1, 2, 3, 4, 5),
  v = c(1, 1, 1, 1, 1)^2
)

N <- 10000
true.mu <- 2
true.sigma_sq <- 3^2

y.obs <- rnorm(N, mean = true.mu, sd = sqrt(true.sigma_sq))

K <- length(parameters$m)
X <- rep(1, K) ## Dirichlet parameter

rigid.hyperprior.data <- list(
     Y = y.obs,
     N = N,
     X = X,
     K = K,
     means = digamma(X)-digamma(X[K]),
     Sigma = constructSigma(X),
     mu = parameters$m,
     sigma_sq = parameters$v)

rigid.dirichlet.posterior <- sampling(compiled.dirichlet, data = rigid.hyperprior.data, refresh = 0,
                               control = list(adapt_delta = .99, max_treedepth = 15))
alphas.rigid.dirichlet <- extract(rigid.dirichlet.posterior, 'alpha')$alpha

print(rigid.dirichlet.posterior, pars = c("mustar", "vstar", "alpha"))
post.alpha.cred.rigid.dirichlet <- apply(alphas.rigid.dirichlet, 2, quantile, probs = c(.025, .975)) 

rigid.logisticNormal.posterior <- sampling(compiled.logisticNormal, data = rigid.hyperprior.data,  refresh = 0,
                                      control = list(adapt_delta = .99, max_treedepth = 15))

print(rigid.logisticNormal.posterior, pars = c("mustar", "vstar", "alpha"))

alphas.rigid.logisticNormal <- extract(rigid.logisticNormal.posterior, 'alpha')$alpha

post.alpha.cred.rigid.logisticNormal <- apply(alphas.rigid.logisticNormal, 2, quantile, probs = c(.025, .975))

KLs <- unlist(lapply(1:K, function(k) {
  kl_gauss(mstar = true.mu, vstar = true.sigma_sq,
           mi = parameters$m[k], vi = parameters$v[k], type = "pf")  
}))
KLs

estimated.alphas <- data.frame(KL_prop = (1/KLs)/sum(1/KLs), dirichlet = colMeans(alphas.rigid.dirichlet), logiNormal = colMeans(alphas.rigid.logisticNormal))
estimated.alphas

pooled.pars <- apply(estimated.alphas, 2, function(x) pool_par_gauss(alpha = x, m = parameters$m, v = parameters$v) )
apply(pooled.pars, 2, function(x) kl_gauss(mstar = true.mu, vstar = true.sigma_sq, mi = x[1], vi = x[2], type = "pf"))

#####

alpha.MC.dirichlet <- rdirichlet(M, X)
alpha.MC.logisticNormal <- rlogisticnorm(N = M,
                                         m = digamma(X)-digamma(X[K]),
                                         Sigma = constructSigma(X))
logkernel <- function(alpha, data){
  y <- data$y
  sigmasq <- data$s2
  ppars <- pool_par_gauss(alpha, m = parameters$m, v = parameters$v)
  m0 <- ppars[1]
  v0 <- ppars[2]
  xbar <- mean(y)
  s2 <-  sum((y-xbar)^2)
  vtilde <- n*v0 + sigmasq 
  l1 <- -s2/(2*sigmasq)
  l2 <- - (n*(xbar-m0)^2)/(2*vtilde)
  l3 <- -0.5 * ( log(vtilde) + (n-1) * log(2*pi*sigmasq) + log(2*pi))
  ans <- l1 + l2 + l3
  return(ans)
}

lZ.dirichlet <- matrixStats::logSumExp(apply(alpha.MC.dirichlet, 1, logkernel, data = list(y = y,  s2 = true.sigma_sq)))-log(M) 
lZ.logisticNormal <- matrixStats::logSumExp(apply(alpha.MC.logisticNormal, 1, logkernel, data = c(y, n)))-log(M) 


phi <- function(alpha, data, lZ){
  y <- data[1]
  n <- data[2]
  ppars <- pool_par(alpha, a = av, b = bv)
  ans <- log(alpha) + lbeta(a = ppars[1] + y, b = ppars[2] + n-y) - lbeta(a = ppars[1], b = ppars[2]) - lZ
  ans <- exp(ans)
  return(ans)
}

round(post.phi.dirichlet <- rowMeans(apply(alpha.MC.dirichlet, 1, phi, data = c(y, n), lZ  = lZ.dirichlet)), 3 )
round(post.phi.logisticNormal <- rowMeans(apply(alpha.MC.logisticNormal, 1, phi, data = c(y, n), lZ = lZ.logisticNormal)), 3)


