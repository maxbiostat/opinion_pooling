source("../pooling_aux.r")
library(matrixStats)
######## Functions
logkernel <- function(alpha, data){
  y <- data$y
  n <- length(y)
  sigmasq <- data$s2
  ppars <- pool_par_gauss(alpha, m = mv, v = vv)
  m0 <- ppars[1]
  v0 <- ppars[2]^2
  ###
  xbar <- mean(y)
  s2 <-  sum((y-xbar)^2)
  vtilde <- n*v0 + sigmasq 
  l1 <- -s2/(2*sigmasq)
  l2 <- - (n*(xbar-m0)^2)/(2*vtilde)
  l3 <- -0.5 * ( log(vtilde) + (n-1) * log(2*pi*sigmasq) + log(2*pi))
  ans <- l1 + l2 + l3
  return(ans)
}
phi <- function(alpha, data, lZ){
  y <- data$y
  n <- length(y)
  sigmasq <- data$s2
  ppars <- pool_par_gauss(alpha, m = mv, v = vv)
  m0 <- ppars[1]
  v0 <- ppars[2]^2
  ###
  xbar <- mean(y)
  s2 <-  sum((y-xbar)^2)
  vtilde <- n*v0 + sigmasq 
  l1 <- -s2/(2*sigmasq)
  l2 <- - (n*(xbar-m0)^2)/(2*vtilde)
  l3 <- -0.5 * ( log(vtilde) + (n-1) * log(2*pi*sigmasq) + log(2*pi))
  ans <- log(alpha) + l1 + l2 + l3 - lZ
  ans <- exp(ans)
  return(ans)
}
get_alphas <- function(mv, vv, y, fixedvar, M, alpha.samples.dirichlet, alpha.samples.logisticNormal){
  K <- length(mv)
  ###$ 1. Compute the  Sample from the prior
  lZ.dirichlet <- matrixStats::logSumExp(apply(alpha.samples.dirichlet, 1,
                                               logkernel, data = list(y = y, s2 = fixedvar)))-log(M) 
  lZ.logisticNormal <- matrixStats::logSumExp(apply(alpha.samples.logisticNormal, 1,
                                                    logkernel, data = list(y = y, s2 = fixedvar)))-log(M) 
  ###$ 2. Compute expectations
  post.phi.dirichlet <- rowMeans(apply(alpha.samples.dirichlet, 1, phi,
                                       data = list(y = y, s2 = sigmasq), lZ  = lZ.dirichlet)) 
  post.phi.logisticNormal <- rowMeans(apply(alpha.samples.logisticNormal, 1, phi,
                                            data = list(y = y, s2 = sigmasq), lZ = lZ.logisticNormal))
  ###$ 3. Compute marginal likelihoods
  log.marginal.likelihoods <- rep(NA, K)
  for (k in 1:K){ 
    log.marginal.likelihoods[k] <- normal_mean_marg_like(y = y, sigma = sqrt(fixedvar),
                                                         m = mv[k], v = vv[k],
                                                         log = TRUE)
  }
  logS <- matrixStats::logSumExp(log.marginal.likelihoods)
  alphapp <- exp(log.marginal.likelihoods-logS)
  return(list(
    dirichlet = post.phi.dirichlet,
    logisticNormal = post.phi.logisticNormal,
    marglike = alphapp
  ))
}
colLwr <- function(x, alpha = .95) apply(x, 2, function(col) quantile(col, probs = (1-alpha)/2 ))
colUpr <- function(x, alpha = .95) apply(x, 2, function(col) quantile(col, probs = (1 +alpha)/2 ))
colMed <- function(x, alpha = .95) apply(x, 2, median)
colSD <- function(x, alpha = .95) apply(x, 2, sd)
compute_posterior_mean_alpha0 <- function(mv, vv, ## priors hyperparameters
                                          truemu, fixedvar, sample_size, ## DGP parameters 
                                          X, ## hyperparameters 
                                          M = 1E4,
                                          n_rep = 100){
  K <- length(mv)
  if(length(vv) != K) stop("Dimensions of mv and vv do not match")
  if(length(X) != K) stop("Dimensions of X and mv do not match")
  ###$ 0.1 Generate data
  ys <- matrix(rnorm(sample_size*nrep, mean = truemu, sd = sqrt(fixedvar)), ncol = n_rep, nrow = sample_size)
  ###$ 0.2 Sample from the prior
  alpha.MC.dirichlet <- LearnBayes::rdirichlet(M, X)
  alpha.MC.logisticNormal <- rlogisticnorm(N = M,
                                           m = digamma(X)-digamma(X[K]),
                                           Sigma = constructSigma(X))
  ###$ 3. run nrep simulations
  simus <- lapply(1:n_rep, function(i){
    get_alphas(mv = mv, vv = vv, y = ys[, i], fixedvar = fixedvar, M = M,
               alpha.samples.dirichlet = alpha.MC.dirichlet,
               alpha.samples.logisticNormal = alpha.MC.logisticNormal)
  })
  
  alphas.dir <- do.call(rbind, lapply(simus, function(x) x$dirichlet))
  alphas.LN <- do.call(rbind, lapply(simus, function(x) x$logisticNormal))
  alphas.MaL  <- do.call(rbind, lapply(simus, function(x) x$marglike))
  ###$ 4. Format output
  out <- data.frame(prior = rep(c("Dirichlet", "Logistic-Normal", "MarginalLikelihoods"), each = K),
                    expert = rep(paste0("expert_", 0:(K-1)), 3),
                    alpha_mean = c(colMeans(alphas.dir),  colMeans(alphas.LN), colMeans(alphas.MaL)),
                    alpha_median = c(colMed(alphas.dir),  colMed(alphas.LN), colMed(alphas.MaL)),
                    alpha_sd = c(colSD(alphas.dir),  colSD(alphas.LN), colSD(alphas.MaL)),
                    alpha_lwr = c(colLwr(alphas.dir),  colLwr(alphas.LN), colLwr(alphas.MaL)),
                    alpha_upr = c(colUpr(alphas.dir),  colUpr(alphas.LN), colUpr(alphas.MaL)),
                    n_experts = K, 
                    n_obs = sample_size,
                    mu = truemu,
                    tau = fixedvar,
                    nrep = n_rep,
                    delta_m = abs(m0-truemu))
}

########
m0 <- 1 ; v0 <- .25^2
m1 <- 2 ; v1 <- .5^2
mv <- c(m0, m1)
vv <- c(v0, v1)
truemu <- 2
sigmasq <- 1^2
X <- rep(2, length(mv))
nrep <- 100

ns <- c(5, 50, 100, 500, 1000, 5000, 10000)
# ns <- c(seq(5, 500, length.out = 10), 1E4)
system.time(
  results.list <- lapply(ns, function(n)
    compute_posterior_mean_alpha0(mv = mv, vv = vv, truemu = truemu, fixedvar = sigmasq, sample_size = n, X = X) 
  )
)

results.dt <- do.call(rbind, results.list)

library(ggplot2)

p0 <- ggplot(results.dt, aes(x = n_obs, y = alpha_mean, col = expert, fill = expert))+
  geom_ribbon(aes(ymin = alpha_lwr, ymax = alpha_upr), alpha = .4) + 
  geom_line() +
  facet_grid(prior~.) +
  theme_bw(base_size = 16)

p1 <- ggplot(subset(results.dt, expert == "expert_0"), aes(x = n_obs, y = alpha_mean,
                                                                 colour = prior, fill = prior))+
  geom_line() +
  scale_x_log10("Sample size (n)", expand = c(0, 0)) +
  scale_y_continuous(expression(alpha[0]), expand = c(0, 0)) +
  geom_ribbon(aes(ymin = alpha_lwr, ymax = alpha_upr), alpha = .4) + 
  theme_bw(base_size = 16)

p2 <- ggplot(subset(results.dt, expert == "expert_0" & prior != "MarginalLikelihoods"),
       aes(x = n_obs, y = alpha_mean, colour = prior, fill = prior))+
  geom_line() +
  scale_x_log10("Sample size (n)", expand = c(0, 0)) +
  scale_y_continuous(expression(alpha[0]), expand = c(0, 0)) +
  geom_ribbon(aes(ymin = alpha_lwr, ymax = alpha_upr), alpha = .4) + 
  theme_bw(base_size = 16)

p1
p2

figW <- 16
figH <- 8

ggsave(filename = "concentration_rates_normal.pdf",
       width = figW, height = figH, units = "cm",
       plot = p1)