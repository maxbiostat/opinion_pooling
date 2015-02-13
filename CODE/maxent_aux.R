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
  c( a[1]/a[2], qgamma( c((1-alpha)/2, (1+alpha)/2), a[1], a[2]))
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
######################################
# Transform + min KL section. These functions only apply to the R_0 example. It's possible to write more general ones, though.
# Computing the induced distribution using pool-then-induce (PI)
# q_1(y) = M(\pi(x) =  dgamma.ratio(y, a1 = sum(alpha*a1p), b1 = sum(alpha*b1p), a2 = sum(alpha*a2p), b2 = sum(alpha*b2p))
# g_i = M(f_i(x)) =  dgamma.ratio(y, k1 = ga1, t1 = gb1, k2 = ga2, t2 = gb2)
########
kl.transform <- function(alpha, a1v, b1v, a2v, b2v, Nv,  ga1, gb1, ga2, gb2, gN){
  kl2int <- function(y){
    dgamma.ratio(y, k1 = ga1, t1 = gb1, k2 = ga2, t2 = gb2, N = gN)* log(
      dgamma.ratio(y, k1 = ga1, t1 = gb1, k2 = ga2, t2 = gb2, N = gN)/
        dgamma.ratio(y, k1 = sum(alpha*a1v), t1 = sum(alpha*b1v),
                     k2 = sum(alpha*a2v), t2 = sum(alpha*b2v), N = sum(alpha*Nv))  
      )
  }  
  integrate(kl2int, 0, Inf)$value
}
