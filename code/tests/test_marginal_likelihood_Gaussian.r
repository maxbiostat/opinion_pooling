### mu ~ normal(m0, v0)
### y|mu ~ normal(mu, tau)
### => y ~ normal(m0, v0 + tau)
truemean <- 2
truevar <- 1^2
m0 <- 0
v0 <- 1^2
n <- 100
# set.seed(666)
y <- rnorm(n = n , mean = truemean, sd = sqrt(truevar))
hist(y, probability = TRUE)

normal.list <- list(
  N = n,
  Y = y,
  tau = truevar,
  m = m0,
  v = v0
)
library(rstan)
normal_model <- stan_model("normal_knownVar.stan")

mcmc <- sampling(normal_model, data = normal.list)

### Approach 0: bridge sampling
MaL.bridge <- bridgesampling::bridge_sampler(mcmc)$logml

### Approach 1: direct calculation

source("../pooling_aux.r")
MaL.direct <- normal_mean_marg_like(y = y, sigma = truevar, m = m0, v = v0, log = TRUE)

MaL.bridge - MaL.direct




