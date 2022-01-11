source("Savchuk_data.R")
library(logPoolR)
#####
get_max_entropy_weights <- function(av, bv, Ntries = 1000, ncores = 2){
  K <- length(av)
  ent.many.startingPoints <- matrix(
    rnorm(n = (K-1)*Ntries, mean = 0, sd = 100),
    ncol = K-1, nrow = Ntries)
  many.ents <- parallel::mclapply(1:Ntries, function(i) {
    optim(ent.many.startingPoints[i, ],
          logPoolR:::optentbeta_inv, # notice this is -H_pi
          ap = av, bp = bv)
  }, mc.cores = ncores)
  optimised.ents <- unlist(lapply(many.ents, function(x) x$value))
  alphaMaxEnt.opt <- logPoolR:::alpha_01(
    many.ents[[which.min(optimised.ents)]]$par
  )
  return(alphaMaxEnt.opt)
}

get_min_KL_weights <- function(av, bv, Ntries = 1000, ncores = 2){
  K <- length(av)
  kl.many.startingPoints <- matrix(rnorm(n = (K-1)*Ntries,
                                         mean = 0, sd = 100),
                                   ncol = K-1, nrow = Ntries)
  many.kls <- parallel::mclapply(1:Ntries, function(i) {
    optim(kl.many.startingPoints[i, ],
          logPoolR:::optklbeta_inv,
          ap = av, bp = bv, type = "fp")
  }, mc.cores = ncores)
  optimised.kls <- unlist(lapply(many.kls, function(x) x$value))
  alphaMinKL.opt <- alpha_01(many.kls[[which.min(optimised.kls)]]$par)
  return(alphaMinKL.opt)
}
### Prep
prior.av <- c(a0, a1, a2, a3) 
prior.bv <- c(b0, b1, b2, b3)
posterior.av <- prior.av + y
posterior.bv <- prior.bv + (n - y)
K <- length(prior.av)
expert.names <- paste0("Expert_", 0:(K-1)) 

### Entropies 
prior.entropies <- posterior.entropies <- rep(NA, K)
for(k in 1:K){
  prior.entropies[k] <- logPoolR:::entropy_beta(prior.av[k], prior.bv[k])
  posterior.entropies[k] <- logPoolR:::entropy_beta(posterior.av[k], posterior.bv[k])
} 

entropy.table <- tibble::tibble(
  expert = expert.names,
  prior_entropy = prior.entropies,
  posterior_entropy = posterior.entropies
)

entropy.table

### Weights

# Equal weights
alpha.equal <- rep(1/K, K)

prior.ab.star.equal <- logPoolR:::pool_par(alpha.equal,
                                           prior.av,
                                           prior.bv)

posterior.ab.star.equal <- logPoolR:::pool_par(alpha.equal,
                                               posterior.av,
                                               posterior.bv)

# BMA
marginal.likelihoods <- rep(NA, K)
for (k in 1:K){
  marginal.likelihoods[k] <- logPoolR:::ml_beta(yi = y,
                                                ni = n,
                                                a = prior.av[k],
                                                b = prior.bv[k])
}

marginal.likelihoods
bma.weights <- marginal.likelihoods/sum(marginal.likelihoods)

posterior.ab.star.bma <- logPoolR::pool_par(alpha = bma.weights,
                                            a = posterior.av,
                                            b = posterior.bv)

# Maximum Entropy
alpha.maxent.prior <- get_max_entropy_weights(av = prior.av,
                                              bv = prior.bv,
                                              ncores = 8)
alpha.maxent.posterior <- get_max_entropy_weights(av = posterior.av,
                                                  bv = posterior.bv,
                                                  ncores = 8)

prior.ab.star.maxent <- pool_par(alpha = alpha.maxent.prior,
                                 a = prior.av,
                                 b = prior.bv)

posterior.ab.star.maxent <- pool_par(alpha = alpha.maxent.posterior,
                                     a = posterior.av,
                                     b = posterior.bv)

# Minimum KL
alpha.minKL.prior <- get_min_KL_weights(av = prior.av,
                                        bv = prior.bv,
                                        ncores = 8)
alpha.minKL.posterior <- get_min_KL_weights(av = posterior.av,
                                            bv = posterior.bv,
                                            ncores = 8)

prior.ab.star.minKL <- pool_par(alpha = alpha.minKL.prior,
                                a = prior.av,
                                b = prior.bv)

posterior.ab.star.minKL <- pool_par(alpha = alpha.minKL.posterior,
                                    a = posterior.av,
                                    b = posterior.bv)

### Hierarchical priors
require("LearnBayes")

M <- 100000
X <- c(1, 1, 1, 1)/10
alpha.MC.Dirichlet <- rdirichlet(M, X)
alpha.MC.logisticNormal <- logPoolR::rlogisticnorm(N = M,
                                                   m = digamma(X)-digamma(X[K]),
                                                   Sigma = logPoolR::constructSigma(X))

apply(alpha.MC.Dirichlet, 2, mean)
apply(alpha.MC.logisticNormal, 2, mean)


prior.ab.star.Dirichlet <- rowMeans(
  apply(alpha.MC.Dirichlet, 1,  function(x){
    logPoolR::pool_par(alpha = x, a = prior.av, b = prior.bv)
  } )
)

apply(alpha.MC.Dirichlet, 2, sd)
apply(alpha.MC.logisticNormal, 2, sd)

prior.ab.star.logisticNormal <- rowMeans(
  apply(alpha.MC.logisticNormal, 1,  function(x){
    logPoolR::pool_par(alpha = x, a = prior.av, b = prior.bv)
  } )
)

### Hierarchical posteriors

library(cmdstanr)
library(rstan)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())


## Dirichlet
compiled.Dirichlet <- cmdstanr::cmdstan_model(
  "../stan/posterior_beta_Dirichlet_pooled.stan")


betadata.stan <- list(Y = y, N = n, K = K,
                      a = prior.av, b = prior.bv,
                      X = X)

Dirichlet.posterior.raw <- compiled.Dirichlet$sample(
  data = betadata.stan,
  chains = 4, parallel_chains = 4,
  adapt_delta = .99,
  max_treedepth = 15,
)
Dirichlet.posterior <- stanfit(Dirichlet.posterior.raw)

print(Dirichlet.posterior, pars = c("theta", "alpha"))

check_hmc_diagnostics(Dirichlet.posterior)

theta.Dirichlet <- extract(Dirichlet.posterior, 'theta')$theta
alphas.Dirichlet <- extract(Dirichlet.posterior, 'alpha')$alpha
post.alpha.cred.Dirichlet <- apply(alphas.Dirichlet, 2,
                                   quantile, probs = c(.025, .975))
betaPars.Dirichlet <- extract(Dirichlet.posterior,
                              c('astar', 'bstar'))
posterior.ab.star.Dirichlet <- unlist(
  lapply(betaPars.Dirichlet, mean)
) 


## Logistic normal
compiled.logisticNormal <- cmdstanr::cmdstan_model(
  "../stan/posterior_beta_logisticNormal_pooled.stan")

betadata.stan$means <- digamma(X)-digamma(X[K])
betadata.stan$Sigma <- constructSigma(X)

logisticNormal.posterior.raw <- compiled.logisticNormal$sample(
  data = betadata.stan,
  chains = 4, parallel_chains = 4,
  adapt_delta = .99,
  max_treedepth = 15,
)
logisticNormal.posterior <- stanfit(logisticNormal.posterior.raw)

print(logisticNormal.posterior, pars = c("theta", "alpha"))

check_hmc_diagnostics(logisticNormal.posterior)

theta.logisticNormal <- extract(logisticNormal.posterior, 'theta')$theta
alphas.logisticNormal <- extract(logisticNormal.posterior, 'alpha')$alpha
post.alpha.cred.logisticNormal <- apply(alphas.logisticNormal, 2,
                                        quantile, probs = c(.025, .975))

betaPars.logisticNormal <- extract(logisticNormal.posterior,
                                   c('astar', 'bstar'))

posterior.ab.star.logisticNormal <- unlist(
  lapply(betaPars.logisticNormal, mean)
)

## Rufo et al. 2012

alphas.rufo <- c(0, 0, 0, 1)

prior.ab.star.Rufo <- logPoolR::pool_par(alpha = alphas.rufo,
                                         a = prior.av,
                                         b = prior.bv)

#### Collating and organizing stuff 

prior.pars <- list(
  equal_weights = prior.ab.star.equal,
  maximum_entropy = prior.ab.star.maxent,
  minimum_KL = prior.ab.star.minKL,
  hierarchical_Dirichlet = prior.ab.star.Dirichlet ,
  hierarchical_LogisticNormal = prior.ab.star.logisticNormal,
  Rufo_2012 = prior.ab.star.Rufo)

lapply(prior.pars,
       function(p){
         logPoolR:::ml_beta(yi = y, ni = n, a = p[1], b = p[2])
       } )

posterior.pars <- list(
  equal_weights = posterior.ab.star.equal,
  maximum_entropy = posterior.ab.star.maxent,
  minimum_KL = posterior.ab.star.minKL,
  hierarchical_Dirichlet = posterior.ab.star.Dirichlet ,
  hierarchical_LogisticNormal = posterior.ab.star.logisticNormal)

save(prior.pars,
     posterior.pars,
     file = "saved_data/hyperparameters_Savchuk.RData")
