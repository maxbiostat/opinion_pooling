pool.par <- function(alpha, a, b){
  # takes an alpha vector and vectors with the K a and b parameters
  # returns a vector with the parameters of the pooled beta/gamma prior
  c(crossprod(a, alpha), crossprod(b, alpha))
}
### Log-pooled density
dpool <- function(x, D, alphas){
  # D is a set (list type) of K of densities
  # alphas is a K-dimensional vector with the weights a_i>0  and sum(alphas)=1
  # if(sum(alphas)!=1){ stop ("The vector of weigths doesn't sum up to unity")}
  if(any(alphas<0)) { stop ("All weigths should be non-negative")}
  sapply(x, function(x) prod(unlist(lapply(D, function(d) d(x)))^alphas))
}
## the constant t(\alpha)
# proper integration
# tpool <- function(alpha, D){
#   1/integrate(function(x){
#     sapply(x, function(e) dpool(e, D = D, alpha = alpha))},
#   -Inf, Inf, rel.tol = .Machine$double.eps^0.3)$value
# }
tpool <- function(alpha, D, lwr = -1E5, upr =  1E5){
  ## shitty trapezoid integration
  1/caTools::trapz(x = seq(lwr, upr, length.out = 1000L),
                   y = dpool(x = seq(lwr, upr, length.out = 1000L), D = D, alpha = alphas))
}
##
tpool.positive <- function(alpha, D){ # For those with positive support
  1/integrate(function(x) sapply(x, function(e) dpool(e, D = D, alpha = alpha)), 0, Inf)$value
} 
#
tpool.unit <- function(alpha, D){ # For those with support in [0, 1]
  1/integrate(function(x) sapply(x, function(e) dpool(e, D = D, alpha = alpha)), 0, .999999999)$value
}
##
makedf <- function(mc.samples){
  kde <- density(mc.samples)
  return(Vectorize(approxfun(kde)))
}
# makedf <- function(mc.samples){
#   kde <- KernSmooth::bkde(mc.samples)
#   return(Vectorize( approxfun(kde)) )
# }
##
dpoolnorm <- function(x, D, alphas){
  # D is a set of K of densities
  # alphas is a K-dimensional vector with the weights a_i>0  and sum(alphas)=1
  dens <- dpool(x, D, alphas)
  return(tpool(alphas, D)*dens)
}
#
dpoolnorm.positive <- function(x, D, alphas){
  # D is a set of K of densities
  # alphas is a K-dimensional vector with the weights a_i>0  and sum(alphas)=1
  dens <- dpool(x, D, alphas)
  return(tpool.positive(alphas, D)*dens)
}
dpoolnorm.unit <- function(x, D, alphas){
  ## same as before, but integrating on the open interval (0, 1)
  dens <- dpool(x, D, alphas)
  return(tpool.unit(alphas, D)*dens)
}
#
dpool.entropy <- function(D, alphas){
  expectlog <- function(x) {-log(x) * dpoolnorm.positive(x, D, alphas)} # using dpoolnorm.positive DELIBERATELY introduces a bug
  ent <- integrate(expectlog, 0, Inf)$value
  return(ent)
}
#
stat.beta <- function(a, alpha = .95){
  # returns the mean and 100*alpha% quantiles
  c( a[1] / (a[1] + a[2]), qbeta( c((1-alpha)/2, (1+alpha)/2), a[1], a[2]))
}
stat.gamma <- function(a, alpha = .95){
  # returns the mean and 100*alpha% quantiles
  c(a[1]/a[2], qgamma( c((1-alpha)/2, (1+alpha)/2), a[1], a[2]))
}
#
fbeta <- function(x, par){
  # returns the density for x \in (0, 1) from a 2-dimensional vector of parameters 
  # par c(a, b) 
  dbeta(x, par[1], par[2])
}
#
fgamma <- function(x, par){
  # returns the density for x \in (0, 1) from a 2-dimensional vector of parameters 
  # par c(a, b) 
  dgamma(x, par[1], par[2])
}
#
entropy.gamma <- function(a, b){
  a - log(b) + log(gamma(a)) + (1-a)*digamma(a)
}
#
entropy.beta <- function(a, b){
  log(beta(a, b)) - (a-1)*digamma(a) - (b-1)*digamma(b) + (a + b -2)*digamma(a+b)
}
#
entropy.normal <- function(v){
  .5*log(2*pi*exp(1)*v)
}
#
ent.surface.beta <- function(av, bv, N = 100){
  # Since the maximum value for astar is max(av) and
  # the same goes for bstar, we can draw a surface for a given set of a's and b's 
  # to look at the face of the entropy surface
  amax <- max(av)
  amin <- min(av)
  bmax <- max(bv)
  bmin <- min(bv)
  as <- seq(amin, amax, length.out = N)
  bs <- seq(bmin, bmax, length.out = N)
  grid <- expand.grid(as, bs)
  Es <- apply(grid, 1, function(row) entropy.beta(row[1], row[2]))
  ME <- matrix(Es, nrow = N)
  return(list(M = ME, as = as, bs = bs))
} 
#
ent.surface.gamma <- function(av, bv, N = 100){
  amax <- max(av)
  amin <- min(av)
  bmax <- max(bv)
  bmin <- min(bv)
  as <- seq(amin, amax, length.out = N)
  bs <- seq(bmin, bmax, length.out = N)
  grid <- expand.grid(as, bs)
  Es <- apply(grid, 1, function(row) entropy.gamma(row[1], row[2]))
  ME <- matrix(Es, ncol = N)
  return(list(M = ME, as = as, bs = bs))
} 
#
gensimplex <- function(K, step = .01, bound = TRUE){
  ## Creates a K-simplex on [0,1]
  ## step is the step size (duh!)
  ## bound is whether to include the boundaries of the simplex (0s and 1s)
  if(bound){ s0 <- seq(0, 1, by = step) } else {s0 <- seq(step, 1-step, by = step)} 
  alphas <- vector(K, mode = "list")
  for(i in 1:K)  alphas[[i]] <- s0
  grid <- expand.grid(alphas)
  pos.ind <- apply(grid, 1, sum)
  simplex <- grid[pos.ind==1, ]
  names(simplex) <- paste("alpha_", 0:(K-1), sep = "")
  return(simplex)
}
#
renorm <- function(x) x + abs(min(x))
#####################################3
## Kullback-Liebler divergence section
kl.gamma <- function(a0, b0, a1, b1){
  ## KL(f||g), where f~G(a0, b0) and g ~ G(a1, b1)
  ## KL(a0, b0, a1, b1)
  (a0-a1)*digamma(a0) - log(gamma(a0)) + log(gamma(a1)) +
    a1*(log(b0/b1)) + a0*((b1-b0)/b0)
}
kl.beta <- function(a0, b0, a1, b1){
  ## KL(f||g), where f ~ B(a0, b0) and g ~ B(a1, b1)
  log(beta(a1, b1)/beta(a0, b0)) + (a0-a1)*digamma(a0) +
    (b0-b1)*digamma(b0) + (a1-a0 + b1-b0)*digamma(a0+b0)
}
#####################################
## optim section
alpha.real <- function(alpha){
  p <- length(alpha)
  log( alpha[-p] / alpha[p]  )
}
##
# alpha.01 <- function(alpha.inv){
#   alpha.inv <- sapply(alpha.inv, function(x) min(x, 10)) ## trying to keep away from boundaries
#   alphap <- 1/(1 + sum(exp(alpha.inv)))
#   alphas <- c(exp(alpha.inv) * alphap, alphap)
#   cat(alpha.inv, "transf: ", alphas, "sum: ", sum(alphas), "\n")
#   if(any(alphas == 1)) alphas[which(alphas != 1)] <- 0 # trying to go round numerical issues
#   return(alphas)
# }
invlogit <- function (x){
  1/(1 + exp(-x))
}
alpha.01 <- function(alpha.inv){
  ## using the representation in the Stan manual
  ## also here https://cran.r-project.org/web/packages/SALTSampler/vignettes/SALTSampler.html
  K <- length(alpha.inv) + 1
  z <- rep(NA, K-1)
  alphas <- rep(0, K)
  for(k in 1:(K-1)){
    z[k] <- invlogit(alpha.inv[k] + log( 1/(K-k) ))
    alphas[k] <-  (1 - sum(alphas[1:(k-1)])) * z[k]
  }
  alphas[K] <- 1-sum(alphas[-K])
  return(alphas)
}
################################
rgelman <- function(N, m, K = NULL, c){
  ## N = number of draws
  ## m vector of means. If a scalar is given, it is recycled K times
  ## c = coefficient of variation. Used to bypass the need to specify a vector of variances. Note this is not in the original paper.
  if(is.null(K)) K <- length(m)
  if(K==1){ means <- rep(m, K)} else{means <- m }  
  draws <- sapply(1:N, function(i) rnorm(K, m = means, s = abs(c*means) ))
  ndraws <- apply(draws, 2, function(x) exp(x)/sum(exp(x)) ) # normalised draws
  return(t(ndraws))
}
################################
constructSigma <- function(X){
  ## Construct the variance-covariance matrix for the logistic normal to match the parameter vector
  ## X, from the Dirichlet distribution
  K <- length(X)
  Sigma <- matrix(rep(trigamma(X[K]), K^2), ncol = K, nrow = K)
  diag(Sigma) <- trigamma(X) + trigamma(X[K])
  return(Sigma)
}
################################
## Let's sample from the prior in Gelman (1995) JCompGraph Stats (http://www.stat.columbia.edu/~gelman/research/published/moments.pdf)
## This is the first example, in page 46. See also Gelman, Bois & Jiang (1996) JASA.
rlogisticnorm <- function(N, m, Sigma){
  ## N = number of draws
  ## m is a vector of means. If a scalar is given, it is recycled K times
  ## Sigma is the variance-covariance matrix
  K <- length(m)
  draws <-  mvtnorm::rmvnorm(n = N, mean = m, sigma = Sigma)
  logisticNorm <- function(x){
    D <- length(x)
    y <- rep(NA, D)
    for(d in 1:(D-1)) y[d] <- exp(x[d])/(1 + sum(exp(x[-D])) )
    y[D] <- 1/(1 + sum(exp(x[-D])) )
    return(y)
  } 
  ndraws <- apply(draws, 1, logisticNorm) # normalised draws
  return(t(ndraws))
}
######
# Logarithmic pooling (LP) via Sampling-importance-resampling (SpIR)
LP.SpIR <- function(k, l, Model, rq1, dq2, dL1, dL2, alpha){
  ## currently tuned to perform LP on two dimensions
  theta.samp <- rq1(n = k)
  phi.transf <- unlist(
    parallel::mclapply(1:nrow(theta.samp), function(i) Model(theta.samp[i, ]), mc.cores = 4) 
  )
  q1star <- makedf(phi.transf)
  getWeight <- function(theta, phi){
    (dq2(phi)/q1star(phi))^{1-alpha} * dL1(theta) * dL2(phi)
  }
  ws <- sapply(seq_len(nrow(theta.samp)), function(i) {
    getWeight(theta = theta.samp[i, ], phi = phi.transf[i])
  })
  ws[which(ws == -Inf)] <- 0 ## giving  weight = 0 to "impossible" values
  resamp.Pos <-  sample(seq_len(nrow(theta.samp)), size = l,
                        replace = TRUE, prob = ws/sum(ws))
  return(list(
    theta.resamp = theta.samp[resamp.Pos, ],
    phi.resamp = phi.transf[resamp.Pos])
  )
}
########################
## grabbed from https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
proposalfunction <- function(param, m, S, h){
  npar <- length(param)
  K <- length(h)
  simplexProp <- SALTSampler::PropStep(y = param[(npar-K + 1):npar], i = sample(1:K, 1), h = h)
  prop <- c(
    mvtnorm::rmvnorm(1, mean = m, sigma = S),
    as.numeric(simplexProp)
  )
  return(
    list(
      value = prop,
      dens.curr = attr(simplexProp, "dbt") * 2, ## should achieve correct Hastings ratio
      dens.prop = attr(simplexProp, "dbt")
    )
  )
}
#
run_metropolis_MCMC <- function(startvalue, iterations,
                                mm, SS, XX, hh,
                                every = 20000, verbose = TRUE){
  chain <- array(dim = c(iterations + 1, length(startvalue)))
  chain[1, ] <- startvalue
  for (i in 1:iterations){
    if(verbose && i%%every==0) print(paste('iteration', i,'of', iterations))
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
##########
LPMCMC <- function(Nit, burnin, startvalue, mm, SS, XX, hh){
  chain <- run_metropolis_MCMC(startvalue = startvalue,
                               mm = mm, SS = SS, X = XX, h = hh,
                               iterations = Nit)
  return(
    list(
      trace = chain[-(1:burnin), ],
      acceptance = 1-mean(duplicated(chain[-(1:burnin), ]))
    )
  )
}