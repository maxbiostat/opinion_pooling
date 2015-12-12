# library(plot3D)
pool.par <- function(alpha, a, b){
  # takes an alpha vector and vectors with the K a and b parameters
  # returns a vector with the parameters of the pooled beta/gamma prior
  c(crossprod(a, alpha), crossprod(b, alpha))
}
### Log-pooled density
dpool <- function(x, D, alphas){
  # D is a set (list type) of K of densities
  # alphas is a K-dimensional vector with the weights a_i>0  and sum(alphas)=1
  #   if(sum(alphas)!=1){ break ("The vector of weigths doesn't sum up to unity")}
  #   if(any(alphas<=0)) { break ("All weigths should be greater than zero")}
  sapply(x, function(x) prod(unlist(lapply(D, function(d) d(x)))^alphas))
}
## the constant t(\alpha)
tpool <- function(alpha, D){
  1/integrate(function(x) sapply(x, function(e) dpool(e, D = D, alpha = alpha)), -Inf, Inf)$value
}
tpool.positive <- function(alpha, D){ # For those with positive support
  1/integrate(function(x) sapply(x, function(e) dpool(e, D = D, alpha = alpha)), 0, Inf)$value
} 
#
tpool.unit <- function(alpha, D){ # For those with positive support
  1/integrate(function(x) sapply(x, function(e) dpool(e, D = D, alpha = alpha)), 0, .999999999)$value
}
##
dpoolnorm <- function(x, D, alphas){
  # D is a set of K of densities
  # alphas is a K-dimensional vector with the weights a_i>0  and sum(alphas)=1
#   if(sum(alphas)!=1){ break ("The vector of weigths doesn't sum up to unity")}
#   if(any(alphas<=0)) { break ("All weigths should be greater than zero")}
  dens <- dpool(x, D, alphas)
  return(tpool(alphas, D)*dens)
}
#
dpoolnorm.positive <- function(x, D, alphas){
  # D is a set of K of densities
  # alphas is a K-dimensional vector with the weights a_i>0  and sum(alphas)=1
#   if(sum(alphas)!=1){ break ("The vector of weigths doesn't sum up to unity")}
#   if(any(alphas<=0)) { break ("All weigths should be greater than zero")}
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
#   if(sum(alphas)!=1){ break ("The vector of weigths doesn't sum up to unity")}
#   if(any(alphas<=0)) { break ("All weigths should be greater than zero")}
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
ent.surface.beta <- function(av, bv){
  # Since the maximum value for astar is max(av) and
  # the same goes for bstar, we can draw a surface for a given set of a's and b's 
  # to look at the face of the entropy surface
  amax <- max(av)
  amin <- min(av)
  bmax <- max(bv)
  bmin <- min(bv)
  as <- seq(amin, amax, length.out = 100)
  bs <- seq(bmin, bmax, length.out = 100)
  grid <- expand.grid(as, bs)
  N <- length(as)
  Es <- apply(grid, 1, function(row) entropy.beta(row[1], row[2]))
  ME <- matrix(Es, nrow = N)
  return(list(M = ME, as = as, bs = bs))
} 
#
ent.surface.gamma <- function(av, bv){
  amax <- max(av)
  amin <- min(av)
  bmax <- max(bv)
  bmin <- min(bv)
  as <- seq(amin, amax, length.out = 100)
  bs <- seq(bmin, bmax, length.out = 100)
  grid <- expand.grid(as, bs)
  N <- length(as)
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
  simplex <- grid[pos.ind==1,]
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

alpha.01 <- function(alpha.inv){
  alphap <- 1/(1 + sum(exp(alpha.inv)))
  c(exp(alpha.inv) * alphap, alphap)
}

################################
## Let's sample from the prior in Gelman (1995) JCompGraph Stats (http://www.stat.columbia.edu/~gelman/research/published/moments.pdf)
## This is the first example, in page 46. See also Gelman, Bois & Jiang (1996) JASA.
rgelman <- function(N, m, K = NULL, c){
  ## N = number of draws
  ## m vector of means. If a scalar is given, it is recycled K times
  ## c = coefficient of variation. Used to bypass the need to specify a vector of variances. Note this is not in the original paper.
  if(is.null(K)) K <- length(m)
  if(length(m)==1){ means <- rep(m, K)}else{means <- m }  
  draws <- sapply(1:N, function(i) rnorm(K, m = means, s = abs(c*means) ))
  ndraws <- apply(draws, 2, function(x) exp(x)/sum(exp(x)) ) # normalised draws
  return(t(ndraws))
}
