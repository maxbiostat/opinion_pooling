#####################################################################################
# f1 is a half-normal with mean 2 and sd .5
# f2 is a half-normal with mean 2.5 and sd 2
# f3 is a Gamma with shape = 16 rate = 8 
# f4 is a log-normal with meanlog = 0.6689426 and sdlog = 0.7033465  
#####################################################################################
f1 <- function(x, log = FALSE){
  ans <- dnorm(x, mean = 2, sd = 0.5, log = TRUE) - pnorm(q = 0, mean = 2, sd = 0.5, log.p = TRUE, lower.tail = FALSE)
  if(!log) ans <- exp(ans)
  return(ans)
}
f1 <- Vectorize(f1)
#
f2 <- function(x, log = FALSE){
  ans <- dnorm(x, mean = 2.5, sd = 2, log = TRUE) - pnorm(q = 0, mean = 2.5, sd = 2, log.p = TRUE, lower.tail = FALSE)
  if(!log) ans <- exp(ans)
  return(ans)
}
f2 <- Vectorize(f2)
#
f3 <- function(x, log = FALSE){
  ans <- dgamma(x, shape = 16, rate = 8, log = log)  
  return(ans)
}
f3 <- Vectorize(f3)
#
f4 <- function(x, log = FALSE){
  ans <- dlnorm(x, meanlog = 0.668942,  sdlog = 0.7033465, log = log)
  return(ans)
}
f4 <- Vectorize(f4)

## Generate alpha
set.seed(215654)
w <- .01
rk <- (1-w)/3
pos <- 2
alpha_raw <- rep(rk, 4)
alpha_raw[pos] <- w
# alpha_raw <- runif(4, 1, 100)
alphas <- alpha_raw/sum(alpha_raw)
alphas

## Create the pooled density, pi(theta)
source('pooling_aux.r')
Ds <- list(f1, f2, f3, f4) # list with the densities
pi_dens_unnorm <- function(x, log = FALSE) dpool(x, D = Ds, alpha = alphas, log = log)
pi_dens <- function(x) dpoolnorm(x, D = Ds, alpha = alphas)

## Samples from individual "posteriors"
M <-  10000
library(truncnorm)
X1 <- rtruncnorm(n = M, mean = 2, sd = .5, a = 0, b = Inf)
X2 <- rtruncnorm(n = M, mean = 2, sd = 2, a = 0, b = Inf)
X3 <- rgamma(n = M, shape = 16, rate = 8)
X4 <- rlnorm(n = M, meanlog =  0.6689426, sdlog = 0.7033465)

## Now we shall sample from the pool directly using Stan
library(rstan)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
pool_mod <- stan_model("stan/sample_from_pool.stan")
simu.data <- list(K = 4, alpha = alphas)
direct.mcmc <- sampling(pool_mod, simu.data, iter = M)
direct.samples <- extract(direct.mcmc, 'theta')$theta


### Trying to approximate the normalising constant 

## stupid method 1, univariate only

hat.f1 <- makedf(mc.samples = X1)
hat.f2 <- makedf(mc.samples = X2)
hat.f3 <- makedf(mc.samples = X3)
hat.f4 <- makedf(mc.samples = X4)

hat.Ds <- list(hat.f1, hat.f2, hat.f3, hat.f4)

approx.unnorm.pi <- function(x){
  ans <- dpool(x, D = hat.Ds, alphas = alphas)
  if(is.na(ans)) ans <- 0
  return(ans)
} 
approx.unnorm.pi <- Vectorize(approx.unnorm.pi)


approx.norm.const <- integrate(approx.unnorm.pi, 0, Inf)$value
log(approx.norm.const)

approx.norm.pi <- function(x, c = approx.norm.const){
  dens <- approx.unnorm.pi(x)
  if(is.na(dens)) dens <- 0
  return(dens/c)
} 
approx.norm.pi  <- Vectorize(approx.norm.pi)

## (Slightly less) stupid method 2: naive importance sampling

functional <- function(x) 1

simple_is <- function(samples, target_dens, prop_dens, parallel = TRUE, cores = 4){
  weight <- function(x) functional(x) * exp( target_dens(x, log = TRUE) - prop_dens(x, log = TRUE))
  weight <- Vectorize(weight)
  if(parallel){
    ests <- unlist(parallel::mclapply(samples, weight, mc.cores = cores) )
  }else{
    ests <- sapply(samples, weight)
  }
  
  return(list(estimates = ests, mu = mean(ests)))
}

system.time(
  IS.ests <- list(
    mu1 = simple_is(samples = X1, target_dens = pi_dens_unnorm, prop_dens = f1),
    mu2 = simple_is(samples = X2, target_dens = pi_dens_unnorm, prop_dens = f2),
    mu3 = simple_is(samples = X3, target_dens = pi_dens_unnorm, prop_dens = f3),
    mu4 = simple_is(samples = X4, target_dens = pi_dens_unnorm, prop_dens = f4)
  )
)

lapply(IS.ests, function(x) log(x$mu))

IS.vars <- unlist(lapply(IS.ests, function(x) var(x$estimates)))
mu.weights <- (1/IS.vars)/sum(1/IS.vars)

## Normalising constant estimates

true.const <- -log(tpool_positive(alpha = alphas, D = Ds)) ## 'true' stuff
bridge.const <- bridgesampling::bridge_sampler(direct.mcmc, silent = TRUE) ## brigdesampling
IS.const <- mean(unlist(lapply(IS.ests, function(x) log(x$mu))))  ## importance sampling
IS.const.wa <- sum(alphas * unlist(lapply(IS.ests, function(x) log(x$mu)))) ## importance sampling (weighted by alpha)
IS.const.wvar <- sum(mu.weights * unlist(lapply(IS.ests, function(x) log(x$mu)))) ## importance sampling (weighted by inv_var)


true.const-bridge.const$logml
true.const-IS.const
true.const-IS.const.wa
true.const-IS.const.wvar

approx.norm.pi.IS <- function(x, c = IS.const, log = FALSE){
  ans <- pi_dens_unnorm(x, log = TRUE) - IS.const.wvar
  if(!log) ans <- exp(ans)
  return(ans)
} 
approx.norm.pi.IS <- Vectorize(approx.norm.pi.IS)

#### The plot

mm <- max(max(X1), max(X2), max(X3), max(X4))

hist(direct.samples, probability = TRUE, xlab = expression(theta), ylim = c(0, .8), breaks = 30, main = "Grey = direct samples")
curve(f1, 0, mm, lwd = 3, add = TRUE)
curve(f2, 0, mm, lwd = 3, col = 2, add = TRUE)
curve(f3, 0, mm, lwd = 3, col = 3, add = TRUE)
curve(f4, 0, mm, lwd = 3, col = 4, add = TRUE)
curve(pi_dens, 0, mm, lwd = 4, lty = 2, col = 5, add = TRUE)
curve(approx.norm.pi, 0, mm, lwd = 4, lty = 3, col = 2, add = TRUE)
curve(approx.norm.pi.IS, 0, mm, lwd = 4, lty = 4, col = 6, add = TRUE)
legend(x = "topright",
       legend = c("f1-normal(2, 0.5)", "f2-normal(2.5, 2)",
                  "f3-gamma(16, 8)", "f4-lognormal(0.67,0.70)", "exact pool", "kernel approx pool", "IS approx pool"),
       col = c(1:5, 2, 6),
       lwd = 2,
       lty = c(rep(1, 4), 2, 3, 4),
       bty = 'n')