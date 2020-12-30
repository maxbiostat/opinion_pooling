source("pooling_aux.r")

m0 <- 1 ; v0 <- 1^2
m1 <- 2 ; v1 <- 10^2
m2 <- 3 ; v2 <- 1^2
m3 <- 4 ; v3 <- 1^2

mv <- c(m0, m1, m2, m3)
vv <- c(v0, v1, v2, v3)
K <- length(mv)

## Observed data
n <- 100
sigmasq <- 2
truemu <- 2
y <- rnorm(n, mean = truemu, sd  = sqrt(sigmasq))

############

library(ggplot2)

mu.grid <- seq(-3*truemu, 3*truemu, length.out = 1000)
expert.densities <- vector(K, mode = "list")
for(k in 1:K){
  expert.densities[[k]] <- data.frame(mu = mu.grid,
                                      dens = dnorm(mu.grid, mean = mv[k], sd = sqrt(vv[k])),
                                      expert = paste("expert_", k-1, sep = ""))
  
}
expert.densities.df <- do.call(rbind, expert.densities)

expert_priors <- ggplot(expert.densities.df, aes(x = mu, y = dens,
                                                 linetype = expert, colour = expert)) + 
  geom_line(size = 2) +
  scale_linetype_manual(values = c("twodash", "dotted", "longdash", "solid"))+
  scale_colour_brewer(palette = "Spectral") +
  scale_x_continuous(expression(mu), expand = c(0, 0)) +
  scale_y_continuous(expression(f[i](mu)), expand = c(0, 0)) +
  theme_bw(base_size = 20)

expert_priors
ggsave(expert_priors, filename = "../plots/expert_densities_SimuNormals.pdf")


####### Hierarchical priors
require("LearnBayes")

M <- 10000
X <- c(1, 1, 1, 1)/10
alpha.MC.dirichlet <- rdirichlet(M, X)
alpha.MC.logisticNormal <- rlogisticnorm(N = M,
                              m = digamma(X)-digamma(X[K]),
                              Sigma = constructSigma(X))

apply(alpha.MC.dirichlet, 2, mean)
apply(alpha.MC.logisticNormal, 2, mean)

apply(alpha.MC.dirichlet, 2, sd)
apply(alpha.MC.logisticNormal, 2, sd)

####### Hierarchical posteriors

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

normaldata.stan <-  list(Y = y, tau = sqrt(sigmasq),
                       X = X, N = n, K = K,
                       mu = mv, sigma_sq = vv)

## Dirichlet
compiled.dirichlet <- stan_model("stan/posterior_fixedVariance_normal_Dirichlet_pooled.stan")
dirichlet.posterior <- sampling(compiled.dirichlet, data = normaldata.stan, 
                                # iter = 5000,
                                control = list(adapt_delta = .99, max_treedepth = 15))
print(dirichlet.posterior, pars = c("theta", "alpha"))

check_hmc_diagnostics(dirichlet.posterior)

mu.dirichlet <- extract(dirichlet.posterior, 'theta')$theta
alphas.dirichlet <- extract(dirichlet.posterior, 'alpha')$alpha
post.alpha.cred.dirichlet <- apply(alphas.dirichlet, 2, quantile, probs = c(.025, .975))
normalPars.dirichlet <- extract(dirichlet.posterior, c('mustar', 'vstar'))
mv.dirichlet <- unlist( lapply(normalPars.dirichlet, mean) ) 
  
## Logistic normal
compiled.logisticNormal <- stan_model("stan/posterior_fixedVariance_normal_logisticNormal_pooled.stan")

normaldata.stan$means <- digamma(X)-digamma(X[K])
normaldata.stan$Sigma <- constructSigma(X)

logisticNormal.posterior <- sampling(compiled.logisticNormal, data = normaldata.stan,
                                     # iter = 5000,
                                     control = list(adapt_delta = .99, max_treedepth = 15))
print(logisticNormal.posterior, pars = c("theta", "alpha"))

check_hmc_diagnostics(logisticNormal.posterior)

mu.logisticNormal <- extract(logisticNormal.posterior, 'theta')$theta
alphas.logisticNormal <- extract(logisticNormal.posterior, 'alpha')$alpha
post.alpha.cred.logisticNormal <- apply(alphas.logisticNormal, 2, quantile, probs = c(.025, .975))
normalPars.logisticNormal <- extract(logisticNormal.posterior, c('mustar', 'vstar'))
mv.logisticNormal <- unlist( lapply(normalPars.logisticNormal, mean) )


#### Testing a simple Monte Carlo approximation of the marginal posterior mean of alpha
logkernel <- function(alpha, data){
  y <- data$y
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

lZ.dirichlet <- matrixStats::logSumExp(apply(alpha.MC.dirichlet, 1, logkernel, data = list(y = y, s2 = sigmasq)))-log(M) 
lZ.logisticNormal <- matrixStats::logSumExp(apply(alpha.MC.logisticNormal, 1, logkernel, data = list(y = y, s2 = sigmasq)))-log(M) 

lZ.dirichlet
bridgesampling::bridge_sampler(dirichlet.posterior)

lZ.logisticNormal
bridgesampling::bridge_sampler(logisticNormal.posterior)

phi <- function(alpha, data, lZ){
  y <- data$y
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

(post.phi.dirichlet <- rowMeans(apply(alpha.MC.dirichlet, 1, phi, data = list(y = y, s2 = sigmasq), lZ  = lZ.dirichlet)) )
(post.phi.logisticNormal <- rowMeans(apply(alpha.MC.logisticNormal, 1, phi, data = list(y = y, s2 = sigmasq), lZ = lZ.logisticNormal)) )

colMeans(alphas.dirichlet)
colMeans(alphas.logisticNormal)
