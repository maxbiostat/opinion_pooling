# Leo Bastos & Luiz Max Carvalho (2017)
source("pooling_aux.r")
source("beta_elicitator.r")

## The idea is as follows: all expert priors but one,  will be rather inconsistent
## with the data.
## Why, you ask? To see how much we can "learn" about \alpha

p0 <- elicit.beta(m0 = 9/10, v0 = 1/100); a0 <- p0$a ; b0 <- p0$b
p1 <- elicit.beta(m0 = 6/10, v0 = 1/100); a1 <- p1$a ; b1 <- p1$b
p2 <- elicit.beta(m0 = 5/10, v0 = 1/100); a2 <- p2$a ; b2 <- p2$b
p3 <- elicit.beta(m0 = 8/10, v0 = 1/100); a3 <- p3$a ; b3 <- p3$b

av <- c(a0, a1, a2, a3)
bv <- c(b0, b1, b2, b3)
K <- length(av)

beta.mode <- function(a, b) (a-1)/(a +b -2)
beta.mean <- function(a, b) a/ (a + b)
beta.sd <- function(a, b) sqrt( (a*b)/( (a + b + 1) * (a + b)^2 ))
beta.mean(av, bv)
beta.sd(av, bv)

# Data
y <- 6
n <- 10
#
alphas <- rep(1/K, K)
###########################
###########################
# Now let's approximate the posterior of theta using a different method
## Priors
beta.PooledPrior <- function(x, weights){
  pstar <- pool.par(alpha = weights, a = av, b = bv)
  dbeta(x, pstar[1], pstar[2], log = TRUE)
} 
#
prior <- function(param, X){
  npar <- length(param) 
  K <- length(X)
  aalphas <- arm::invlogit(param[(npar - K + 1):npar])
  pr <- beta.PooledPrior(arm::invlogit(param[1]), weights = aalphas) + 
    log(gtools::ddirichlet(aalphas, alpha = X)) 
  return(pr)
}
## Likelihoods
Lik <- function(param){
  dbinom(x = y, size = n, prob = arm::invlogit(param[1]), log = TRUE)
}
## Posterior
posterior <- function(param, X, verbose = TRUE){
  Lk <- Lik(param)
  pr <-  prior(param, X = X)
  if(verbose) cat("Likelihood:", Lk, " prior:", pr,
                  " theta:", arm::invlogit(param[1]), 
                  " alphas:", alpha.01(param[2:K]),
                  "\n")
  return (Lk + pr)
}
#
## Compile stuff to try and improve performance
prior <- compiler::cmpfun(prior)
Lik <- compiler::cmpfun(Lik)
posterior <- compiler::cmpfun(posterior)
##
Niter <- 5000
run_metropolis_MCMC <- function(startvalue, Niter, X){
  require(MHadaptive)
  chain <- Metro_Hastings(li_func = posterior, X = X, pars = startvalue,
                          iterations = Niter, burn_in = round(.2 * Niter),
                          par_names = c('theta', paste('alpha', 0:(K-1))),
                          quiet = FALSE)
  return(chain)
}
system.time(
  Chain <- run_metropolis_MCMC(
    startvalue = arm::logit(c(.5, alphas)),
    X = rep(1/2, K),
    Niter = Niter)
)
Chain$acceptance
#
MH.Post <- data.frame(theta = arm::invlogit(Chain$trace[, 1]),
                      alphas = data.frame(t(apply(Chain$trace[, -1], 1, alpha.01))))
#
apply(MH.Post, 2, coda::effectiveSize)
apply(MH.Post, 2,  function(x) c(mean = mean(x), quantile(x, probs = c(.025, .5, .975))))

norml <- function(x) x/sum(x)
shouldBe.alpha <- norml(1/(abs(beta.mean(av, bv)- y/n) + .1)) ## closer to data gets more weight
post.pstar <- pool.par(shouldBe.alpha, av, bv) + c(y, n-y)
hist(MH.Post$theta, probability = TRUE)
curve(dbeta(x, shape1 = post.pstar[1], shape2 = post.pstar[2]), lwd = 2, add = TRUE)