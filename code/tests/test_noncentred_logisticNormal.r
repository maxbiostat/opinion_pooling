source("../pooling_aux.r")
source("../beta_elicitator.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

compiled.logisticNormal <- stan_model("../stan/posterior_beta_logisticNormal_pooled.stan")
compiled.logisticNormalNC <- stan_model("../stan/posterior_beta_logisticNormalNonCentred_pooled.stan")

x <- 5
n <- 10
cv <- .2

pars <- get_parameter_vectors(cv_correct = cv)
K <- length(pars$av)
X <- rep(1, K)/10

betadata.stan <- list(Y = x, X = X,
                      means = digamma(X)-digamma(X[K]),
                      Sigma = constructSigma(X),
                      nu = 1,
                      s = 1,
                      N = n, K = K, a = pars$av, b = pars$bv)



logisticNormal.posterior <- sampling(compiled.logisticNormal, data = betadata.stan,
                                     # refresh = 0,
                                     control = list(adapt_delta = .99))


logisticNormalNC.posterior <- sampling(compiled.logisticNormalNC, data = betadata.stan,
                                     # refresh = 0,
                                     control = list(adapt_delta = .99))

logisticNormal.posterior
logisticNormalNC.posterior
