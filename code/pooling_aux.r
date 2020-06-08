######################## Pooling-related functions #############
pool_par <- function(alpha, a, b){
  if(any(alpha <0)) stop("weights must be non-negative")
  if(!identical(sum(alpha), 1)) stop("weights do not sum to 1")
  # takes an alpha vector and vectors with the K a and b parameters
  # returns a vector with the parameters of the pooled beta/gamma prior
  c(crossprod(a, alpha), crossprod(b, alpha))
}
pool_par_gauss <- function(alpha, m, v){
  ## same as pool_par(), but for the Gaussian distribution
  ## Takes in MEANS and VARIANCES and outputs MEAN and SD (for plotting reasons)
  ws <- alpha/v
  vstar <-  1/sum(ws)
  mstar <- sum(ws*m) * vstar
  c(mstar, sqrt(vstar))
}
### Log-pooled density
dpool <- function(x, D, alphas, log = FALSE){
  # D is a set (list type) of K of densities
  # alphas is a K-dimensional vector with the weights a_i>0  and sum(alphas)=1
  # if(sum(alphas)!=1){ stop ("The vector of weigths doesn't sum up to unity")}
  if(any(alphas<0)) { stop ("All weigths should be non-negative")}
  log_dens <- function(x) sum(unlist(lapply(D, function(d) log(d(x)) )) * alphas)
  ans <- sapply(x,  function(x) sum(log_dens(x)))
  if(!log) ans <- exp(ans)
  return(ans)
}
### the constant t(\alpha)
tpool <- function(alpha, D, trapez = TRUE, lwr = -1E5, upr =  1E5){
  if(trapez){
    the_int <- caTools::trapz(x = seq(lwr, upr, length.out = 1000L),
                              y = dpool(x = seq(lwr, upr, length.out = 1000L), D = D, alphas = alpha))
  }else{
    # proper integration
    the_int <- integrate(function(x){
      sapply(x, function(e) dpool(e, D = D, alphas = alpha))},
      -Inf, Inf)$value 
  }
  return(1/the_int)
}
##
tpool_positive <- function(alpha, D, trapez = FALSE, lwr = 0, upr =  1E5){
  if(trapez){
    ## shitty trapezoid integration
    the_int <- caTools::trapz(x = seq(lwr, upr, length.out = 1000L),
                              y = dpool(x = seq(lwr, upr, length.out = 1000L), D = D, alphas = alpha))
  }else{
    ## proper integration
    the_int <- integrate(function(x) dpool(x, D = D, alphas = alphas), 0, Inf)$value
  }
  return(1/the_int)
}
##
tpool_unit <- function(alpha, D){ # For those with support in [0, 1]
  1/integrate(function(x) sapply(x, function(e) dpool(e, D = D, alphas = alpha)), 0, .999999999)$value
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
  return(tpool_positive(alphas, D)*dens)
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

### General simplex functions

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
renorm <- function(x) x + abs(min(x))

alpha_real <- function(alpha){
  p <- length(alpha)
  log( alpha[-p] / alpha[p]  )
}
invlogit <- function (x){
  1/(1 + exp(-x))
}
alpha_01 <- function(alpha.inv){
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
LP_SpIR <- function(k, l, Model, rq1, dq2, dL1, dL2, alpha, cores = 4){
  require(compiler)
  ## currently tuned to perform LP on two dimensions
  theta.samp <- rq1(n = k)
  phi.transf <- unlist(
    parallel::mclapply(1:nrow(theta.samp), function(i) Model(theta.samp[i, ]), mc.cores = cores) 
  )
  q1star <- makedf(phi.transf)
  get_Weight <- function(theta, phi, log = FALSE){
    # (dq2(phi)/q1star(phi))^{1-alpha} * dL1(theta) * dL2(phi)
    ans <- (1-alpha)*{log(dq2(phi))-log(q1star(phi)) } + log(dL1(theta)) + log(dL2(phi))  
    if(!log) ans <- exp(ans)
    return(ans)
  }
  get_Weight <- cmpfun(get_Weight)
  ws <- unlist(
    parallel::mclapply(seq_len(nrow(theta.samp)), function(i) {
      get_Weight(theta = theta.samp[i, ], phi = phi.transf[i])
    }, mc.cores = cores)
  ) 
  ws[which(ws == -Inf)] <- 0 ## giving  weight = 0 to "impossible" values
  ws[which(ws == "NaN")] <- 0 ## giving  weight = 0 to "weird" values ;-)
  resamp.Pos <-  sample(seq_len(nrow(theta.samp)), size = l,
                        replace = TRUE, prob = ws/sum(ws))
  return(list(
    theta.resamp = theta.samp[resamp.Pos, ],
    phi.resamp = phi.transf[resamp.Pos])
  )
}
##
LPSpIR_varying_alpha <- function(k, l, Model, rq1, dq2, dL1, dL2, cores = 4){
  ## currently tuned to perform LP on two dimensions
  theta.samp <- rq1(n = k)
  phi.transf <- unlist(
    parallel::mclapply(1:nrow(theta.samp), 
                       function(i) Model(theta.samp[i, 1:2]), mc.cores = cores) 
  )
  q1star <- makedf(phi.transf)
  corrected_q1star <-  function(x){
    res <- q1star(x)
    res[is.na(res)] <- 0
    return(res)
  } 
  getKa <- function(a){
    tpool(alpha = c(a, 1-a),
          D = list(
            f0 = function(x) corrected_q1star(x),
            f1 = function(x) dnorm(x, mean = 7800, sd = 1300) )
    )
  }
  getKa <- compiler::cmpfun(getKa)
  #
  get_Weight <- function(theta, phi, log = FALSE){
    # getKa(theta[3]) * (dq2(phi)/q1star(phi))^{1-theta[3]} * dL1(theta) * dL2(phi)
    ans <- log(getKa(theta[3])) + (1-theta[3])*{log(dq2(phi)) - log(q1star(phi))} + log(dL1(theta)) + log(dL2(phi))
    if(!log) ans <- exp(ans)
    return(ans)
  }
  get_Weight <- compiler::cmpfun(get_Weight)
  ws <- unlist(
    parallel::mclapply(seq_len(nrow(theta.samp)), function(i) {
      get_Weight(theta = theta.samp[i, ], phi = phi.transf[i])
    }, mc.cores = cores)
  ) 
  ws[which(ws == -Inf)] <- 0 ## giving  weight = 0 to "impossible" values
  ws[which(ws == "NaN")] <- 0 ## giving  weight = 0 to "weird" values ;-)
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
######################## General utils #############
mean_ci <- function(x, alpha = .95, na.rm = TRUE){
  qs <- quantile(x, probs = c(.5 * (1 - alpha), .5 * (1 + alpha) ), na.rm = na.rm)
  return(data.frame(
    mean = mean(x, na.rm = na.rm),
    lwr = as.numeric(qs[1]),
    upr = as.numeric(qs[2])
  ))
}
get_ratio <- function(x) {
  ## outputs the ratio of the largest and second largest values in a vector x
  y <- x[order(x, decreasing = TRUE)]
  y[1]/y[2]
}
#### Plotting
require(ggplot2)
coord_radar <- function (theta = "x", start = 0, direction = 1) 
{ ## from http://web-r.org/board_ISVF22/8271
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") 
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}
number_ticks <- function(n) {function(limits) pretty(limits, n)}
######################## Distributions #############
### Beta distribution
ml_beta <- function(yi, ni, a, b, log = FALSE){ # Equation (9) in Raftery et al (2007) 
  ans <- lgamma(ni + 1) - lgamma(ni - yi + 1) - lgamma(yi + 1)  +
    lgamma(a + b)-lgamma(a + b + ni) +
    ( lgamma(a+yi)-lgamma(a) ) + ( lgamma(b + ni - yi) -lgamma(b))
  if(!log) ans <- exp(ans)
  return(ans)
}
#
marginal_likelihood_surface_beta <- function(y, n, av, bv, N = 100){
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
  MLs <- apply(grid, 1, function(row) ml_beta(yi = y, ni = n, a = row[1], b = row[2]))
  ML <- matrix(MLs, nrow = N)
  return(list(M = ML, as = as, bs = bs))
} 
#
stat_beta <- function(a, alpha = .95){
  # returns the mean and 100*alpha% quantiles
  c( a[1] / (a[1] + a[2]), qbeta( c((1-alpha)/2, (1+alpha)/2), a[1], a[2]))
}

beta_mode <- function(a, b) (a-1)/(a +b -2)
beta_mean <- function(a, b) a/ (a + b)
beta_sd <- function(a, b) sqrt( (a*b)/( (a + b + 1) * (a + b)^2 ))

fbeta <- function(x, par){
  # returns the density for x \in (0, 1) from a 2-dimensional vector of parameters 
  # par c(a, b) 
  dbeta(x, par[1], par[2])
}
entropy_beta <- function(a, b){
  lbeta(a, b) - (a-1)*digamma(a) - (b-1)*digamma(b) + (a+b-2)*digamma(a+b)
}
optentbeta <- function(alpha, ap, bp){
  entropy_beta(a = sum(alpha*ap), b = sum(alpha*bp))
}
#
optentbeta_inv <- function(alpha.inv, ap, bp){
  alpha <- alpha_01(alpha.inv)
  -optentbeta(alpha, ap, bp)
}
entropy_surface_beta <- function(av, bv, N = 100){
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
  Es <- apply(grid, 1, function(row) entropy_beta(a = row[1], b = row[2]))
  ME <- matrix(Es, nrow = N)
  return(list(M = ME, as = as, bs = bs))
} 
kl_beta <- function(astar, bstar, ai, bi, type = c("pf", "fp")){
  ## if type = pf, computes KL(pi||f),
  ## if type = fp, computes KL(f || pi)
  # where pi ~ Beta(astar, bstar) and f ~ B(ai, bi)
  if(type == "pf"){
    a1 = astar
    b1 = bstar
    a2 = ai
    b2 = bi  
  }else{
    a1 = ai
    b1 = bi
    a2 = astar
    b2 = bstar  
  }
  res <-  lbeta(a2, b2) - lbeta(a1, b1) + (a1-a2)*digamma(a1) +
    (b1-b2)*digamma(b1)  + (a2 - a1 + b2 - b1)*digamma(a1 + b1)
  return(res)
}
optklbeta <- function(alpha, ap, bp, type){
  K <- length(alpha)
  astar <- sum(alpha*ap)
  bstar <- sum(alpha*bp)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){ ds[i] <-  kl_beta(astar = astar, bstar = bstar,
                                    ai = ap[i], bi = bp[i], type = type)} 
  return(ds)
}
optklbeta_inv <- function(alpha.inv, ap, bp, type = "pf"){
  alpha <- alpha_01(alpha.inv)
  sum(optklbeta(alpha, ap, bp, type))
}

### Gamma distribution
stat_gamma <- function(a, alpha = .95){
  # returns the mean and 100*alpha% quantiles
  c(a[1]/a[2], qgamma( c((1-alpha)/2, (1 + alpha)/2), a[1], a[2]))
}
fgamma <- function(x, par){
  # returns the density for x \in (0, 1) from a 2-dimensional vector of parameters 
  # par c(a, b) 
  dgamma(x, par[1], par[2])
}
entropy_gamma <- function(a, b){
  a - log(b) + log(gamma(a)) + (1-a)*digamma(a)
}
entropy_surface_gamma <- function(av, bv, N = 100){
  amax <- max(av)
  amin <- min(av)
  bmax <- max(bv)
  bmin <- min(bv)
  as <- seq(amin, amax, length.out = N)
  bs <- seq(bmin, bmax, length.out = N)
  grid <- expand.grid(as, bs)
  Es <- apply(grid, 1, function(row) entropy_gamma(a = row[1], b = row[2]))
  ME <- matrix(Es, ncol = N)
  return(list(M = ME, as = as, bs = bs))
}
kl_gamma <- function(astar, bstar, ai, bi, type = c("pf", "fp")){
  ## KL(f||g), where f~G(a0, b0) and g ~ G(a1, b1)
  ## KL(a, b0, a1, b1)
  if(type == "pf"){
    a0 = astar
    b0 = bstar
    a1 = ai
    b1 = bi  
  }else{
    a0 = ai
    b0 = bi
    a1 = astar
    b1 = bstar  
  }
  ans <- (a0-a1)*digamma(a0) - log(gamma(a0)) + log(gamma(a1)) +
    a1*(log(b0/b1)) + a0*((b1-b0)/b0)
  return(ans)
}
optklgamma <- function(alpha, ap, bp, type){
  K <- length(alpha)
  astar <- sum(alpha*ap)
  bstar <- sum(alpha*bp)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){ ds[i] <-  kl_gamma(astar = astar, bstar = bstar,
                                     ai = ap[i], bi = bp[i], type = type)} 
  return(ds)
}
optklgamma_inv <- function(alpha.inv, ap, bp, type = "pf"){
  alpha <- alpha_01(alpha.inv)
  sum(optklgamma(alpha, ap, bp, type))
}


### Gaussian distribution
stat_gauss <- function(p, alpha = .95){
  # returns the mean and 100*alpha% quantiles
  c( p[1], qnorm( c((1-alpha)/2, (1+alpha)/2), p[1], p[2]))
}
entropy_gauss <- function(m, v){
  .5*log(2*pi*exp(1)*v)
}
### WARNING: there's no need to optimise the entropy, as there is a "closed-form" solution, namely picking the distribution with 
### the largest variance, but we'll indulge for the sake of completeness.
optentgauss <- function(alpha, mp, vp){
  ws <- alpha/vp
  mstar <- sum(ws*mp)/sum(ws)
  vstar <-  1/sum(ws)
  entropy_gauss(m = mstar, v = vstar)
}
#
optentgauss_inv <- function(alpha.inv, mp, vp){
  alpha <- alpha_01(alpha.inv)
  -optentgauss(alpha, mp, vp)
}
kl_gauss <- function(mstar, vstar, mi, vi, type = c("pf", "fp")){
  ## WARNING: parametrisation is MEAN and VARIANCE!
  if(type == "pf"){
    m0 = mstar
    v0 = vstar
    m1 = mi
    v1 = vi  
  }else{
    m0 = mi
    v0 = vi
    m1 = mstar
    v1 = vstar  
  }
  ans <- log(sqrt(v1)) - log( sqrt(v0))  + (v0 + (m0-m1)^2)/(2 * v1) - 1/2
  return(ans)
}
optklgauss <- function(alpha, mp, vp, type){
  K <- length(alpha)
  ws <- alpha/vp
  mstar <- sum(ws*mp)/sum(ws)
  vstar <-  1/sum(ws)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){ ds[i] <-  kl_gauss(mstar = mstar, vstar = vstar,
                                     mi = mp[i], vi = vp[i], type = type)} 
  return(ds)
}
optklgauss_inv <- function(alpha.inv, mp, vp, type = "pf"){
  alpha <- alpha_01(alpha.inv)
  sum(optklgauss(alpha, mp, vp, type))
}
