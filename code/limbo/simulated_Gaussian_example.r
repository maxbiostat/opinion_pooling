source("pooling_aux.r")
library(rstan)
rstan_options(auto_write = TRUE)

compiled.dirichlet <- stan_model("stan/posterior_normal_Dirichlet_pooled.stan")
compiled.logisticNormal <-   stan_model("stan/posterior_normal_logisticNormal_pooled.stan")


###################
parameters <- list(
  m = c(1, 2, 3, 4, 5),
  v = c(5, 4, 3, 2, 1)^2
)

N <- 1000
true.mu <- 1
true.sigma_sq <- 1^2


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

