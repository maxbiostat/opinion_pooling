# Leo Bastos & Luiz Max Carvalho (2017)
source("pooling_aux.r")
source("beta_elicitator.r")

meta <- read.csv("../data/meta_analysis_Malta_2010.csv")

meta$SampleSize

meta$av <- meta$HIV + 1
meta$bv <- meta$SampleSize - meta$HIV + 1

K <- nrow(meta)

beta.mode(meta$av, meta$bv)
beta.mean(meta$av, meta$bv)
beta.sd(meta$av, meta$bv)


# Entropy surface (for a future dominance analysis) 
library(fields)
ES <- ent.surface.beta(meta$av, meta$bv)
image.plot(ES$as, ES$bs, ES$M, xlab = expression(a), ylab = expression(b))

# log-Pooled prior: Beta(a'alpha, b'alpha)

# Preparing the table 
PaperBeta.tbl <- data.frame( mean.prior = rep(NA, 5), lower.prior = NA, upper.prior = NA)
rownames(PaperBeta.tbl) <- c("Equal weights", "maxEnt", "min KL div.", "Hier. prior Diri", "Hier prior Exp")

# AlphasBeta.tbl <- data.frame(matrix(NA, nrow = 4, ncol = length(meta$av)))
# rownames(AlphasBeta.tbl) <- c("maxEnt", "min KL div.", "Hier. prior Diri", "Hier prior Exp")
AlphasBeta.tbl <- data.frame(matrix(NA, nrow = 2, ncol = length(meta$av)))
rownames(AlphasBeta.tbl) <- c("maxEnt", "min KL div.")
colnames(AlphasBeta.tbl) <- paste("alpha", 0:(K-1), sep = "")

######################################################
###### Equal weights alphas
######################################################
alphaEqual <- rep(1/K, K)

ab.Equal.star <- pool.par(alphaEqual, meta$av, meta$bv)

# Prior
(PaperBeta.tbl[1, 1:3] <- stat.beta(ab.Equal.star))


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

a.ent <- optim(c(0,0,0,0,0), optentbeta.inv, ap = meta$av, bp = meta$bv,
               method = "SANN", control=list(maxit = 100000))
alphaMaxEnt.opt <- alpha.01(a.ent$par)
(AlphasBeta.tbl[1,] <- alphaMaxEnt.opt)

ab.MaxEnt.star <- pool.par(alphaMaxEnt.opt, meta$av, meta$bv)

# Prior
( PaperBeta.tbl[2, 1:3] <- stat.beta(ab.MaxEnt.star) )


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

a.kl <- optim(c(0, 0, 0, 0, 0), optklbeta.inv, ap = meta$av, bp = meta$bv)
#            method = "SANN", control=list(maxit = 100000))

alphaKL.opt <- alpha.01(a.kl$par)
(AlphasBeta.tbl[2, ] <- alphaKL.opt)

ab.KL.star <- pool.par(alphaKL.opt, meta$av, meta$bv)

# Prior
(PaperBeta.tbl[3, 1:3] <- stat.beta(ab.KL.star))


######################################################
###### Hierarchical prior
# \pi(theta|alpha)
# alpha ~ Dirichlet (X)
# X = {x_0, x_1, ..., x_K}
######################################################
require("LearnBayes")

M <- 100000
X <- c(1, 1, 1, 1, 1, 1)
alpha.MC.dir <- rdirichlet(M, X)
alpha.MC.exp <- rlogisticnorm(N = M,
                              m = digamma(X)-digamma(X[K]),
                              Sigma = constructSigma(X))

apply(alpha.MC.dir, 2, mean)
apply(alpha.MC.exp, 2, mean)

apply(alpha.MC.dir, 2, sd)
apply(alpha.MC.exp, 2, sd)

beta.par.dir <- alpha.MC.dir %*% cbind(meta$av, meta$bv)
beta.par.exp <- alpha.MC.exp %*% cbind(meta$av, meta$bv)

theta.par.dir <- apply(beta.par.dir, 1, function(x) rbeta(1, x[1], x[2]))
theta.par.exp <- apply(beta.par.exp, 1, function(x) rbeta(1, x[1], x[2]))

# Prior
PaperBeta.tbl[4, 1] <- mean(theta.par.dir)
PaperBeta.tbl[4, 2:3] <- quantile(theta.par.dir, c(.025, .975))

PaperBeta.tbl[5, 1] <- mean(theta.par.exp)
PaperBeta.tbl[5, 2:3] <- quantile(theta.par.exp, c(.025, .975))


######################################################
############## Results
######################################################
# Table

round(PaperBeta.tbl, 2)
round(AlphasBeta.tbl, 2)
###  Plotting
# png("../manuscript/figures/beta_example_simulated_data.png")
# par(mfrow = c(2, 1))
ccx <- 1.5
curve(fbeta(x, par = c(meta$av[1], meta$bv[1]) ), 0, .5,  ylab = "Density", main = "Study induced distributions",
      xlab = expression(theta), lwd = 3 , lty = 1, cex.lab = ccx, cex.axis = ccx, cex.main = ccx, cex.sub = ccx)
for( k in 2:K)
  curve(fbeta(x, par = c(meta$av[k], meta$bv[k]) ), 0, .5, lwd = 3, col = k, lty = k, add = TRUE)

legend(x = "topright",
       legend = paste("Study", 1:K),
       col = 1:K, lwd = 3, lty = 1:K, bty = "n"
)
# dev.off()
# Combined Priors
curve(fbeta(x, par = ab.Equal.star), 0, 0.5, ylab = "Density", main = "log-Pooled distributions",
      xlab = expression(theta), lwd = 3, lty = 1,
      cex.lab = ccx, cex.axis = ccx, cex.main = ccx, cex.sub = ccx)
curve(fbeta(x, par = ab.MaxEnt.star), 0, 1, col = 2, add = TRUE, lwd = 3, lty = 2)
curve(fbeta(x, par = ab.KL.star), 0, 1, col = 3, add = TRUE, lwd = 3, lty = 3)
lines(density(theta.par.dir), col = 4, lwd = 3, lty = 4)
lines(density(theta.par.exp), col = 5, lwd = 3, lty = 5)
legend(x = "topright", bty = "n", col = 1:5, lty = 1:5, cex = 0.8, lwd = 3,
       legend = c("Equal weights (1/K)", "MaxEnt",
                  "MinKL", "Hierarchical Dirichlet", "Hierarchical logistic-normal")
       )
# pars <- list(equal = ab.Equal.star,
#              max_entropy = ab.MaxEnt.star,
#              KL = ab.KL.star,
#              hierD = ab.Hier.star.dir,
#              hierE = ab.Hier.star.exp)
# lapply(pars, function(p) ml.beta(yi = y, ni = n, a = p[1], b = p[2]))