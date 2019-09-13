# Leo Bastos & Luiz Max Carvalho (2017)
source("pooling_aux.r")
source("beta_elicitator.r")

p0 <- elicit.beta(m0 = 7/10, v0 = 1/100); a0 <- p0$a ; b0 <- p0$b
p1 <- elicit.beta(m0 = 7/10, v0 = 1/10); a1 <- p1$a ; b1 <- p1$b
p2 <- elicit.beta(m0 = 7/10, v0 = 1/5); a2 <- p2$a ; b2 <- p2$b
p3 <- elicit.beta(m0 = 7/10, v0 = 1/50); a3 <- p3$a ; b3 <- p3$b

av <- c(a0, a1, a2, a3)
bv <- c(b0, b1, b2, b3)
K <- length(av)

beta.mode <- function(a, b) (a-1)/(a +b -2)
beta.mean <- function(a, b) a/ (a + b)
beta.sd <- function(a, b) sqrt( (a*b)/( (a + b + 1) * (a + b)^2 ))
beta.mean(av, bv)
beta.sd(av, bv)
# Entropy surface (for a future dominance analysis) 
library(fields)
ES <- ent.surface.beta(av, bv)
image.plot(ES$as, ES$bs, ES$M, xlab = expression(a), ylab = expression(b))

# Data
y <- 70
n <- 100

### Marginal (integrated) likelihoods for each of the experts' priors

ml.beta <- function(yi, ni, a, b){ # Equation (9) in Raftery et al (2007) 
  ( gamma(ni + 1)/{gamma(ni - yi + 1) * gamma(yi + 1)} ) *
    ( gamma(a + b)/gamma(a + b + ni)  ) *
    ( gamma(a+yi)/gamma(a) ) * ( gamma(b + ni - yi)/gamma(b)) 
}

marglikes <- rep(NA, K)
for (k in 1:K){ marglikes[k] <- ml.beta(yi = y, ni = n, a = av[k], b = bv[k]) }

marglikes

# log-Pooled prior: Beta(a'alpha, b'alpha)

# Preparing the table 
PaperBeta.tbl <- data.frame( mean.prior = rep(NA, 5), lower.prior = NA, 
                             upper.prior = NA, mean.post = NA, lower.post = NA, 
                             upper.post = NA)
rownames(PaperBeta.tbl) <- c("Equal weights", "maxEnt", "min KL div.", "Hier. prior Diri", "Hier prior Exp")

AlphasBeta.tbl <- data.frame(matrix(NA, nrow = 4, ncol = length(av)))
rownames(AlphasBeta.tbl) <- c("maxEnt", "min KL div.", "Hier. prior Diri", "Hier prior Exp")
colnames(AlphasBeta.tbl) <- paste("alpha", 0:(length(av)-1), sep = "")
######################################################
###### Equal weights alphas
######################################################
alphaEqual <- rep(1/K, K)

ab.Equal.star <- pool.par(alphaEqual, av, bv)
# Prior
(PaperBeta.tbl[1, 1:3] <- stat.beta(ab.Equal.star))
# Posterior
(PaperBeta.tbl[1, 4:6] <- stat.beta(ab.Equal.star + c(y, n - y)))

######################################################
##### Maximizing Entropy
# Maximise H(\pi; alpha), the entropy of the pooled prior
######################################################
optentbeta <- function(alpha, ap, bp){
  entropy.beta(a = sum(alpha*ap), b = sum(alpha*bp))
}

optentbeta.inv <- function(alpha.inv, ap, bp){
  alpha <- alpha.01(alpha.inv)
  -optentbeta(alpha, ap, bp)
}

a.ent <- optim(c(0, 0, 0), optentbeta.inv, ap = av, bp = bv) 
#            method = "SANN", control=list(maxit = 100000))
alphaMaxEnt.opt <- alpha.01(a.ent$par)
(AlphasBeta.tbl[1,] <- alphaMaxEnt.opt)

ab.MaxEnt.star <- pool.par(alphaMaxEnt.opt, av, bv)

# Prior
(PaperBeta.tbl[2, 1:3] <- stat.beta(ab.MaxEnt.star))

# Posterior
(PaperBeta.tbl[2, 4:6] <- stat.beta(ab.MaxEnt.star+ c(y, n - y)))

######################################################
# Minimising KL divergence between each density and the pooled prior
# F = { f_0, f_1, ..., f_K}
# d_i = KL(f_i || \pi) ; L(alpha) = sum(d_i)
# minimise the loss function L(F; alpha)
######################################################
optklbeta <- function(alpha, ap, bp){
  K <- length(alpha)
  astar <- sum(alpha*ap)
  bstar <- sum(alpha*bp)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){ ds[i] <-  kl.beta(astar = astar, bstar = bstar, ai = ap[i], bi = bp[i])} 
  return(ds)
}
optklbeta.inv <- function(alpha.inv, ap, bp){
  alpha <- alpha.01(alpha.inv)
  sum(optklbeta(alpha, ap, bp))
}

a.kl <- optim(c(0, 0, 0), optklbeta.inv, ap = av, bp = bv)
#            method = "SANN", control=list(maxit = 100000))

alphaKL.opt <- alpha.01(a.kl$par)
(AlphasBeta.tbl[2, ] <- alphaKL.opt)

ab.KL.star <- pool.par(alphaKL.opt, av, bv)

# Prior
(PaperBeta.tbl[3, 1:3] <- stat.beta(ab.KL.star))

# Posterior
(PaperBeta.tbl[3, 4:6] <- stat.beta(ab.KL.star + c(y, n-y)))

######################################################
###### Hierarchical prior
# \pi(theta|alpha)
# alpha ~ Dirichlet (X)
# X = {x_0, x_1, ..., x_K}
######################################################
require("LearnBayes")

M <- 10000
X <- c(1, 1, 1, 1)/2 # Jeffreys' prior
alpha.MC.dir <- rdirichlet(M, X)
alpha.MC.exp <- rlogisticnorm(N = M,
                              m = digamma(X)-digamma(X[K]),
                              Sigma = constructSigma(X))

apply(alpha.MC.dir, 2, mean)
apply(alpha.MC.exp, 2, mean)

apply(alpha.MC.dir, 2, sd)
apply(alpha.MC.exp, 2, sd)

apply(alpha.MC.dir, 2, quantile, probs = .975)
apply(alpha.MC.exp, 2, quantile, probs = .975)


beta.par.dir <- alpha.MC.dir %*% cbind(av, bv)
beta.par.exp <- alpha.MC.exp %*% cbind(av, bv)

theta.par.dir <- apply(beta.par.dir, 1, function(x) rbeta(1, x[1], x[2]))
theta.par.exp <- apply(beta.par.exp, 1, function(x) rbeta(1, x[1], x[2]))
# Prior
PaperBeta.tbl[4, 1] <- mean(theta.par.dir)
PaperBeta.tbl[4, 2:3] <- quantile(theta.par.dir, probs = c(.025, .975))

PaperBeta.tbl[5, 1] <- mean(theta.par.exp)
PaperBeta.tbl[5, 2:3] <- quantile(theta.par.exp, probs = c(.025, .975))

# Posterior 
library(rstan)
betadata.stan <-  list(Y = y, X = X, N = n, K = K, a = av, b = bv)
betadata.stan.exp <-  list(Y = y,
                           means = digamma(X)-digamma(X[K]),
                           Sigma = constructSigma(X),
                           N = n, K = K, a = av, b = bv)
# hierpost.dir <- stan(file = "posterior_beta_Dirichlet_pooled.stan",
#                      data = betadata.stan, iter = 1, thin = 1, chains = 1)
# save(hierpost.dir, file = "compiled_beta_post_Dirichlet_sampler.RData")
# hierpost.exp <- stan(file = "posterior_beta_logisticNormal_pooled.stan",
#                      data = betadata.stan.exp, iter = 1, thin = 1, chains = 1)
# save(hierpost.exp, file = "compiled_beta_post_logisticNormal_sampler.RData")
load("compiled_beta_post_Dirichlet_sampler.RData")
load("compiled_beta_post_logisticNormal_sampler.RData")
hierpostsamp.dir <- stan(fit= hierpost.dir,
                     data = betadata.stan, iter = 5000, thin = 1, chains = 1)
hierpostsamp.exp <- stan(fit= hierpost.exp,
                         data = betadata.stan.exp, iter = 5000, thin = 1, chains = 1)

posteriors.dir <- extract(hierpostsamp.dir)
posteriors.exp <- extract(hierpostsamp.exp)
alphas.exp <- matrix(NA, nrow = nrow(posteriors.exp$m), ncol = K )
for (i in 1: nrow(posteriors.exp$m)){ alphas.exp[i, ] <- exp(posteriors.exp$m[i, ])/ sum(exp(posteriors.exp$m[i, ])) }

PaperBeta.tbl[4, 4] <- mean(posteriors.dir$theta)
PaperBeta.tbl[4, 5:6] <- quantile(posteriors.dir$theta, c(.025, .975))

PaperBeta.tbl[5, 4] <- mean(posteriors.exp$theta)
PaperBeta.tbl[5, 5:6] <- quantile(posteriors.exp$theta, c(.025, .975))

post.alpha.dir <- apply(posteriors.dir$alpha, 2, mean)
post.alpha.exp <- apply(alphas.exp, 2, mean)

post.alpha.dir.cred <- apply(posteriors.dir$alpha, 2, quantile, probs = c(.025, .975))
post.alpha.exp.cred <- apply(alphas.exp, 2, quantile, probs = c(.025, .975))

(AlphasBeta.tbl[3, ] <- post.alpha.dir)
(AlphasBeta.tbl[4, ] <- post.alpha.exp)

ab.Hier.star.dir <- pool.par(post.alpha.dir, av, bv)
ab.Hier.star.exp <- pool.par(post.alpha.exp, av, bv)

######################################################
############## Results
######################################################
###  Plotting
# png("../manuscript/figures/beta_example_simulated_data.png")
# par(mfrow = c(2, 1))
ccx <- 1.5
curve(fbeta(x, par = c(a0, b0) ), .5, 1,  ylab = "Density", main = "Expert Priors",
      xlab = expression(theta), lwd = 3 , lty = 3, cex.lab = ccx, cex.axis = ccx, cex.main = ccx, cex.sub = ccx)
curve(fbeta(x, par = c(a1, b1) ), .5, 1, lwd = 3, col = 1, lty = 4, add = TRUE)
curve(fbeta(x, par = c(a2, b2) ), .5, 1, lwd = 3, col = 1, lty = 5, add = TRUE)
curve(fbeta(x, par = c(a3, b3) ), .5, 1, lwd = 3, col = 1, lty = 6, add = TRUE)
legend(x = "topleft",
       legend = paste("Expert", 1:4),
       col = 1, lwd = 3, lty = 3:6, bty = "n"
)
# dev.off()
# Combined Priors
curve(fbeta(x, par = ab.Equal.star), .5, 1, ylab = "Density", main = "Pooled Priors and Posteriors",
      xlab = expression(theta), lwd = 2, lty = 2,
      cex.lab = ccx, cex.axis = ccx, cex.main = ccx, cex.sub = ccx)
curve(fbeta(x, par = ab.MaxEnt.star), .5, 1, col = 2, add = TRUE, lwd = 2, lty = 2)
curve(fbeta(x, par = ab.KL.star), .5, 1, col = 3, add = TRUE, lwd = 2, lty = 2)
lines(density(theta.par.dir), col = 4, lwd = 2, lty = 2)
lines(density(theta.par.exp), col = 5, lwd = 2, lty = 2)
# png("../manuscript/figures/beta_example_simulated_data.png")
# Posteriors
curve(fbeta(x, par = ab.Equal.star+ c(y, n-y)), .5, 1, ylab = "Density", main = "",
      xlab = expression(theta), lwd = 2, add = TRUE)
curve(fbeta(x, par = ab.MaxEnt.star+ c(y, n-y)), .5, 1, col = 2, add = TRUE, lwd = 2)
curve(fbeta(x, par = ab.KL.star+ c(y, n-y)), .5, 1, col = 3, add = TRUE, lwd = 2)
lines(density(posteriors.dir$theta), xlim = c(.5, 1), col = 4, lwd = 2)
lines(density(posteriors.exp$theta), xlim = c(.5, 1), col = 5, lwd = 2)
legend(x = "topleft", bty = "n", col = 1:5, cex = .7,
       legend = c("Equal weights (1/K)", "MaxEnt",
                  "MinKL", "Hierarchical Dirichlet", "Hierarchical ExpNormal"),
       lwd = 2, lty = 1)
legend(x = "top", bty = "n", col = 1,
       legend = c("Priors", "Posteriors"),
       lwd = 2, lty = 2:1)
abline( v = y/n, lwd = 3, lty = 3)
# dev.off()
#############
Posterior_experts <- data.frame(
  alpha = c(as.numeric(c(AlphasBeta.tbl[1, ], AlphasBeta.tbl[2, ])),
            post.alpha.dir, post.alpha.exp),
  lwr = c(rep(NA, 8), post.alpha.dir.cred[1, ], post.alpha.exp.cred[1, ]),
  upr = c(rep(NA, 8), post.alpha.dir.cred[2, ], post.alpha.exp.cred[2, ]),
  expert = rep(paste("expert_", 0:(K-1), sep = ""), 2),
  prior = rep(c("maxEnt", "minKL", "Prior_Dirichlet", "Prior_NormalExp"), each = K)
)
library(ggplot2)
##
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
####
ggplot(data = Posterior_experts,
       aes(x = expert, y = alpha, group = prior, colour = prior, fill = prior)) +
  geom_point() +
  geom_polygon(alpha = 0.4) +
  theme_bw() + 
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
#############
# Now  let's look at marginal likelihoods for the pooled priors

pars <- list(equal = ab.Equal.star,
             max_entropy = ab.MaxEnt.star,
             KL = ab.KL.star,
             hierD = ab.Hier.star.dir,
             hierE = ab.Hier.star.exp)
lapply(pars, function(p) ml.beta(yi = y, ni = n, a = p[1], b = p[2]))


round(PaperBeta.tbl, 2)
round(AlphasBeta.tbl, 2)

get_odds <- function(x) {
  y <- x[order(x, decreasing = TRUE)]
  y[1]/y[2]
}
apply(AlphasBeta.tbl, 1, get_odds)