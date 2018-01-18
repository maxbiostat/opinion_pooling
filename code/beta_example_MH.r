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
## Priors
beta.PooledPrior <- function(x, weights){
  pstar <- pool.par(alpha = weights, a = av, b = bv)
  dbeta(x, pstar[1], pstar[2], log = TRUE)
} 
#
invlogitJac <- function(y){
  arm::invlogit(y) * (1-arm::invlogit(y))
}
alpha01Jac <- function(alpha.inv){
  ## using the representation (and formulae) in the Stan manual
  K <- length(alpha.inv) + 1
  z <- part <- rep(NA, K-1)
  alphas <- rep(0, K)
  for(k in 1:(K-1)){
    z[k] <- invlogit(alpha.inv[k] + log( 1/(K-k) ))
    alphas[k] <-  (1 - sum(alphas[1:(k-1)])) * z[k]
    part[k] <- (1-z[k]) * alphas[k]
  }
  return(prod(part))
}
#
prior <- function(param, X){
  npar <- length(param)
  K <- length(X)
  aalphas <- alpha.01(param[(npar-K + 2):npar])
  beta.PooledPrior(invlogit(param[1]), weights = aalphas) + log(invlogitJac(param[1])) +
    log(gtools::ddirichlet(aalphas, alpha = X)) + alpha01Jac(param[(npar-K + 2):npar])
}
## Likelihoods
Lik <- function(param){
  dbinom(x = y, size = n, prob = invlogit(param[1]), log = TRUE) 
}
## Posterior
posterior <- function(param, X, verbose = FALSE){
  Lk <- Lik(param)
  pr <-  prior(param, X = X)
  if(verbose) cat("Likelihood:", Lk, " prior:", pr,
                  " theta:", invlogit(param[1]), 
                  " alphas:", invlogit(param[2:K]),
                  "\n")
  return (Lk + pr)
}
## grabbed from https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
proposalfunction <- function(param, m, S, h){
  npar <- length(param)
  K <- length(h)
  simplexProp <- SALTSampler::PropStep(y = param[(npar - K + 1):npar],
                                       i = sample(1:K, 1), h = h)
  prop <- c(
    rnorm(1, mean = m, sd = S),
    as.numeric(simplexProp)
  ) 
  return(
    list(
      value = prop,
      dens.curr = attr(simplexProp, "dbt"), 
      dens.prop = attr(simplexProp, "dbt")* 2 ## should achieve correct Hastings ratio
    )
  )
}
#
run_metropolis_MCMC <- function(startvalue, iterations,
                                mm, SS, XX, hh,
                                every = 1000, verbose = TRUE){
  chain <- array(dim = c(iterations + 1, length(startvalue)))
  chain[1, ] <- startvalue
  for (i in 1:iterations){
    if(verbose && i%%every==0) print(paste('iteration', i, 'of', iterations))
    proposal <- proposalfunction(chain[i, ], m = mm, S = SS, h = hh)
    probab  <- exp( (posterior(proposal$value, X = XX) +  proposal$dens.prop) -
                      (posterior(chain[i, ], X = XX) + proposal$dens.curr) )
    if (runif(1) < probab){
      chain[i+1, ] <- proposal$value
    }else{
      chain[i+1, ] <- chain[i, ]
    }
  }
  return(chain)
}
##
LPMCMC <- function(Nit, burnin, startvalue, mm, SS, XX, hh){
  chain <- run_metropolis_MCMC(startvalue = startvalue,
                               mm = mm, SS = SS, X = XX, h = hh,
                               iterations = Nit)
  return(
    list(
      trace = chain,
      acceptance = 1-mean(duplicated(chain[-(1:burnin), ]))
    )
  )
}
## Compile stuff to try and improve performance
prior <- compiler::cmpfun(prior)
Lik <- compiler::cmpfun(Lik)
posterior <- compiler::cmpfun(posterior)
#
proposalfunction <- compiler::cmpfun(proposalfunction)
run_metropolis_MCMC <- compiler::cmpfun(run_metropolis_MCMC)
##
Niter <- 50000
system.time(
  Chain <- LPMCMC(startvalue = c(0, SALTSampler::Logit(rep(1/K, K))),
                  mm = 0,
                  SS = 1,
                  XX = rep(1/2, K), ## Jeffreys prior
                  hh = rep(.2, K),
                  Nit = Niter, burnin = round(.5 * Niter))
)
Chain$acceptance

Transform <- function(x){
    invlogit(x)
}
Posterior <- data.frame(matrix(apply(Chain$trace, 1, Transform), ncol = 5, byrow = TRUE))
names(Posterior) <- c("theta", paste("alpha", 0:(K-1), sep = "_"))  
apply(Posterior, 2, coda::effectiveSize)
apply(Posterior, 2,  function(x) c(mean = mean(x), quantile(x, probs = c(.025, .5, .975))))

plot(coda::as.mcmc(Posterior))

norml <- function(x) x/sum(x)
( shouldBe.alpha <- norml(1/(abs(beta.mean(av, bv)- y/n) + 1/n))  ) ## closer to data gets more weight
post.pstar <- pool.par(colMeans(Posterior[, -1]), av, bv) + c(y, n-y)
shouldBe.pstar <- pool.par(shouldBe.alpha, av, bv) + c(y, n-y)
hist(Posterior$theta, probability = TRUE, xlab = expression(theta))
abline(v = y/n, lwd = 3, lty = 2)
curve(dbeta(x, shape1 = post.pstar[1], shape2 = post.pstar[2]), lwd = 2, add = TRUE)
curve(dbeta(x, shape1 = shouldBe.pstar[1], shape2 = shouldBe.pstar[2]),
      lwd = 2, lty = 2, col = 2,  add = TRUE)
legend(x = "topleft", legend = c("ideal alpha", "posterior alpha"), 
       col = 2:1, lwd = 2, bty = 'n')
