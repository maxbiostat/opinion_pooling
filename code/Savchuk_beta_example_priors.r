# This example was taken from Savchuk and Martz (1994)

source("pooling_aux.r")

a0 <- 18.1 ; b0 <- .995
a1 <- 3.44 ; b1 <- .860 
a2 <- 8.32 ; b2 <- .924
a3 <- 1.98 ; b3 <- .848

av <- c(a0, a1, a2, a3)
bv <- c(b0, b1, b2, b3)
K <- length(av)

# Individual entropies
entropies <- rep(NA, K)
for(k in 1:K) entropies[k] <- entropy_beta(av[k], bv[k])

entropies

## Entropy surface (for a future dominance analysis) 

library(fields)
ES <- entropy_surface_beta(av, bv)
export <- TRUE
if(export){
  pdf("../plots/entropy_surface_failureProbExample.pdf")  
}
image.plot(ES$as, ES$bs, ES$M,
           xlab = expression(a), ylab = expression(b), horizontal = TRUE,
           cex.lab = 1.5, cex.axis = 1.5, axis.args = list(font = 2),
           legend.cex = 1.5,
           legend.lab = expression(H[pi]), main = "Entropy Beta distribution", font = 2)
if(export){
  dev.off()  
}

## Observed data
y <- 9
n <- 10

## Marginal likelihoods

marginal.likelihoods <- rep(NA, K)
for (k in 1:K){ marginal.likelihoods[k] <- ml_beta(yi = y, ni = n, a = av[k], b = bv[k]) }

marginal.likelihoods
round( normalised.marginal.likelihoods<- marginal.likelihoods/sum(marginal.likelihoods), 2 )


MLS <- marginal_likelihood_surface_beta(y = y, n = n, av = av, bv = bv)
export <- TRUE
if(export){
  pdf("../plots/marginalLikelihood_surface_failureProbExample.pdf")  
}
image.plot(MLS$as, MLS$bs, MLS$M,
           xlab = expression(a), ylab = expression(b), horizontal = TRUE,
           cex.lab = 1.5, cex.axis = 1.5, axis.args = list(font = 2),
           legend.cex = 1.5,
           legend.lab = expression(l(y, n)), main = "Marginal likelihood Beta distribution", font = 2)
if(export){
  dev.off()  
}

############
PaperBeta.tbl <- data.frame(mean.prior = rep(NA, 6), lower.prior = NA, 
                             upper.prior = NA, mean.post = NA, lower.post = NA, 
                             upper.post = NA)
rownames(PaperBeta.tbl) <- c("equal_weights", "maximum_entropy", "minimum_KL",
                             "hierarchical_Dirichlet", "hierarchical_LogisticNormal", "Rufo_2012")

AlphasBeta.tbl <- data.frame(matrix(NA, nrow = 5, ncol = length(av)))
rownames(AlphasBeta.tbl) <- c("maximum_entropy", "minimum_KL",
                              "hierarchical_Dirichlet", "hierarchical_LogisticNormal", "Rufo_2012")
colnames(AlphasBeta.tbl) <- paste("alpha_", 0:(K-1), sep = "")


library(ggplot2)

theta.grid <- seq(0, 1, length.out = 1000)
expert.densities <- vector(K, mode = "list")
for(k in 1:K){
  expert.densities[[k]] <- data.frame(theta = theta.grid,
                                      dens = dbeta(theta.grid, shape1 = av[k], shape2 = bv[k]),
                                      expert = paste("expert_", k-1, sep = ""))
  
}
expert.densities.df <- do.call(rbind, expert.densities)

expert_priors <- ggplot(expert.densities.df, aes(x = theta, y = dens,
                                                 linetype = expert, colour = expert)) + 
  geom_line(size = 2) +
  scale_linetype_manual(values = c("twodash", "dotted", "longdash", "solid"))+
  scale_colour_brewer(palette = "Spectral") +
  scale_x_continuous(expression(theta), expand = c(0, 0)) +
  scale_y_continuous(expression(f[i](theta)), expand = c(0, 0)) +
  theme_bw(base_size = 20)

expert_priors
ggsave(expert_priors, filename = "../plots/expert_densities_Savchuk.pdf")

###### Equal weights

alphaEqual <- rep(1/K, K)

ab.Equal.star <- pool_par(alphaEqual, av, bv)
# Prior
(PaperBeta.tbl[1, 1:3] <- stat_beta(ab.Equal.star))
# Posterior
(PaperBeta.tbl[1, 4:6] <- stat_beta(ab.Equal.star + c(y, n - y)))

####### Maximum entropy

N <- 1000 ## could increase to, say, 10000 in order to make sure, but it's fine
ent.many.startingPoints <- matrix(rnorm(n = (K-1)*N, mean = 0, sd = 100), ncol = K-1, nrow = N)
many.ents <- lapply(1:N, function(i) {
  optim(ent.many.startingPoints[i, ], optentbeta_inv, ap = av, bp = bv)
})
optimised.ents <- unlist(lapply(many.ents, function(x) x$value))

hist(optimised.ents)
abline(v = optimised.ents[which.min(optimised.ents)], lty = 2, lwd = 2)

alphaMaxEnt.opt <- alpha_01(many.ents[[which.min(optimised.ents)]]$par)

round(alphaMaxEnt.opt, 2)

( AlphasBeta.tbl[1, ] <- alphaMaxEnt.opt )

ab.MaxEnt.star <- pool_par(alphaMaxEnt.opt, av, bv)

# Prior
(PaperBeta.tbl[2, 1:3] <- stat_beta(ab.MaxEnt.star))

# Posterior
(PaperBeta.tbl[2, 4:6] <- stat_beta(ab.MaxEnt.star +  c(y, n - y)))

####### Minimum KL

N <- 1000 ## could increase to, say, 10000 in order to make sure, but it's fine
kl.many.startingPoints <- matrix(rnorm(n = (K-1)*N, mean = 0, sd = 100), ncol = K-1, nrow = N)
many.kls <- lapply(1:N, function(i) {
  optim(kl.many.startingPoints[i, ], optklbeta_inv, ap = av, bp = bv, type = "fp")
})
optimised.kls <- unlist(lapply(many.kls, function(x) x$value))

hist(optimised.kls)
abline(v = optimised.kls[which.min(optimised.kls)], lty = 2, lwd = 2)

alphaKL.opt <- alpha_01(many.kls[[which.min(optimised.kls)]]$par)

(AlphasBeta.tbl[2, ] <- alphaKL.opt)

ab.KL.star <- pool_par(alphaKL.opt, av, bv)

# Prior
(PaperBeta.tbl[3, 1:3] <- stat_beta(ab.KL.star))

# Posterior
(PaperBeta.tbl[3, 4:6] <- stat_beta(ab.KL.star + c(y, n-y)))
  
####### Hierarchical priors
require("LearnBayes")

M <- 100000
X <- c(1, 1, 1, 1)/10
alpha.MC.dirichlet <- rdirichlet(M, X)
alpha.MC.logisticNormal <- rlogisticnorm(N = M,
                              m = digamma(X)-digamma(X[K]),
                              Sigma = constructSigma(X))

apply(alpha.MC.dirichlet, 2, mean)
apply(alpha.MC.logisticNormal, 2, mean)

apply(alpha.MC.dirichlet, 2, sd)
apply(alpha.MC.logisticNormal, 2, sd)

beta.par.dirichlet <- alpha.MC.dirichlet %*% cbind(av, bv)
beta.par.logisticNormal <- alpha.MC.logisticNormal %*% cbind(av, bv)

theta.par.dirichlet <- apply(beta.par.dirichlet, 1, function(x) rbeta(1, x[1], x[2]))
theta.par.logisticNormal <- apply(beta.par.logisticNormal, 1, function(x) rbeta(1, x[1], x[2]))
# Prior
PaperBeta.tbl[4, 1] <- mean(theta.par.dirichlet)
PaperBeta.tbl[4, 2:3] <- quantile(theta.par.dirichlet, c(.025, .975))

PaperBeta.tbl[5, 1] <- mean(theta.par.logisticNormal)
PaperBeta.tbl[5, 2:3] <- quantile(theta.par.logisticNormal, c(.025, .975))

####### Hierarchical posteriors

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

betadata.stan <-  list(Y = y, X = X, N = n, K = K, a = av, b = bv)

## Dirichlet
compiled.dirichlet <- stan_model("stan/posterior_beta_Dirichlet_pooled.stan")
dirichlet.posterior <- sampling(compiled.dirichlet, data = betadata.stan, 
                                control = list(adapt_delta = .99, max_treedepth = 15))
print(dirichlet.posterior, pars = c("theta", "alpha"))

check_hmc_diagnostics(dirichlet.posterior)

theta.dirichlet <- extract(dirichlet.posterior, 'theta')$theta
alphas.dirichlet <- extract(dirichlet.posterior, 'alpha')$alpha
post.alpha.cred.dirichlet <- apply(alphas.dirichlet, 2, quantile, probs = c(.025, .975))
betaPars.dirichlet <- extract(dirichlet.posterior, c('astar', 'bstar'))
ab.dirichlet <- unlist( lapply(betaPars.dirichlet, mean) ) 
  
(PaperBeta.tbl[4, 4:6] <- mean_ci(theta.dirichlet) )
(AlphasBeta.tbl[3, ] <- colMeans(alphas.dirichlet)) 

## Logistic normal
compiled.logisticNormal <- stan_model("stan/posterior_beta_logisticNormal_pooled.stan")

betadata.stan$means <- digamma(X)-digamma(X[K])
betadata.stan$Sigma <- constructSigma(X)

logisticNormal.posterior <- sampling(compiled.logisticNormal, data = betadata.stan,
                                     control = list(adapt_delta = .99, max_treedepth = 15))
print(logisticNormal.posterior, pars = c("theta", "alpha"))

check_hmc_diagnostics(logisticNormal.posterior)

theta.logisticNormal <- extract(logisticNormal.posterior, 'theta')$theta
alphas.logisticNormal <- extract(logisticNormal.posterior, 'alpha')$alpha
post.alpha.cred.logisticNormal <- apply(alphas.logisticNormal, 2, quantile, probs = c(.025, .975))
betaPars.logisticNormal <- extract(logisticNormal.posterior, c('astar', 'bstar'))
ab.logisticNormal <- unlist( lapply(betaPars.logisticNormal, mean) )

( PaperBeta.tbl[5, 4:6] <- mean_ci(theta.logisticNormal) )
( AlphasBeta.tbl[4, ] <- colMeans(alphas.logisticNormal) ) 

####### KL prior from Rufo et al 2012

alphas.rufo <- c(0, 0, 0, 1)

(AlphasBeta.tbl[5, ] <- alphas.rufo)

ab.rufo <- pool_par(alphas.rufo, av, bv)

# Prior
( PaperBeta.tbl[6, 1:3] <- stat_beta(ab.rufo) )

# Posterior
( PaperBeta.tbl[6, 4:6] <- stat_beta(ab.rufo +  c(y, n - y)) )

#### Finally, tables!
round(PaperBeta.tbl, 3)
round(AlphasBeta.tbl, 3)

round(PaperBeta.tbl, 2)
round(AlphasBeta.tbl, 2)
####### Plotting

posterior_experts <- data.frame(
  alpha = as.numeric(c(AlphasBeta.tbl[1, ], AlphasBeta.tbl[2, ], AlphasBeta.tbl[5, ],
            AlphasBeta.tbl[3, ], AlphasBeta.tbl[4, ])),
  lwr = c(rep(NA, 12), post.alpha.cred.dirichlet[1, ], post.alpha.cred.logisticNormal[1, ]),
  upr = c(rep(NA, 12), post.alpha.cred.dirichlet[2, ], post.alpha.cred.logisticNormal[2, ]),
  expert = rep(paste("expert_", 0:(K-1), sep = ""), 5),
  method = rep(c("maximum_entropy", "minimum_KL", "Rufo_2012",  "Dirichlet", "logistic_normal"), each = K)
)

####
radar_alphas <- ggplot(data = posterior_experts,
       aes(x = expert, y = alpha, group = method, colour = method, fill = method)) +
  geom_point() +
  geom_polygon(alpha = 0.4) +
  theme_bw(base_size = 16) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                     breaks = number_ticks(10)) + 
  coord_radar() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )
radar_alphas

ggsave(plot = radar_alphas, filename = "../plots/alphas_radar_Savchuk.pdf")


#############
# Now  let's look at marginal likelihoods for the pooled priors

pars <- list(equal_weights = ab.Equal.star,
             maximum_entropy = ab.MaxEnt.star,
             minimum_KL = ab.KL.star,
             hierarchical_Dirichlet = ab.dirichlet ,
             hierarchical_LogisticNormal = ab.logisticNormal,
             Rufo_2012 = ab.rufo)
lapply(pars, function(p) ml_beta(yi = y, ni = n, a = p[1], b = p[2]))

apply(AlphasBeta.tbl, 1, get_ratio)

J <- length(pars)
method.densities.list <- vector(J, mode ="list")
for (j in 1:J){
  method.densities.list[[j]] <- data.frame(
    theta = theta.grid,
    dens = dbeta(theta.grid, shape1 = pars[[j]][1], shape2 = pars[[j]][2]),
    method = names(pars)[j]
  )
}

method.densities.df <- do.call(rbind, method.densities.list)

method_priors <- ggplot(method.densities.df, aes(x = theta, y = dens,
                                                        linetype = method, colour = method)) + 
  geom_line(size = 2) +
  scale_x_continuous(expression(theta), expand = c(0, 0)) +
  scale_y_continuous(expression(pi(theta)), expand = c(0, 0)) +
  geom_vline(xintercept = y/n, linetype = "dashed", size = 1.2) +
  theme_bw(base_size = 16)

method_priors
ggsave(method_priors, filename = "../plots/method_prior_densities_Savchuk.pdf")