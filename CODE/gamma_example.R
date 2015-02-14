 source("maxent_aux.R")

  ## Five experts give their parameters for gamma distributions
  ## n = 5 observations from a Poisson (lamb) are made
  set.seed(8340279)
 
  a0 <- 2 ; b0 <- 1.1
  a1 <- 5 ; b1 <- 5
  a2 <- 3 ; b2 <- 6
  a3 <- 1.8 ; b3 <- 1.2
  a4 <- 1.6; b4 <- .8
  
  av <- c(a0, a1, a2, a3, a4)
  bv <- c(b0, b1, b2, b3, b4)
  K <- length(av)
  
  # Data
  lamb <- 1 # 'true' lambda 
  n <- 5
  y <- rpois(n, lamb) 
  ys <- sum(y)
 
 Prior.exps <- data.frame(Prior_0 = stat.gamma(c(a0, b0)),
                          Prior_1 = stat.gamma(c(a1, b1)),
                          Prior_2 = stat.gamma(c(a2, b2)), 
                          Prior_3 = stat.gamma(c(a3, b3)), 
                          Prior_4 = stat.gamma(c(a4, b4)))
 rownames(Prior.exps) <- c("Mean", "Lower", "Upper")
  
 ###########
 ml.gamma <- function(y, a, b){
   n <- length(y)
   s <- sum(y)
   num <- ((b + n)^ (a + s) )*gamma(a + s)
   den <-  (b^a) * gamma(a)
   p <- factorial(prod(y))
   return( (num/den) * 1/p)
 }
 
 marglikes <- rep(NA, K)
 for (k in 1:K){ marglikes[k] <- ml.gamma(y = y, a = av[k], b = bv[k]) }
  marglikes
 
 
  # log-Pooled prior: Gamma(a'alpha, b'alpha)
  
  # Preparing the table 
  PaperGamma.tbl <- data.frame(mean.prior = rep(NA, 4), lower.prior = NA, 
                               upper.prior = NA, mean.post = NA, lower.post = NA, 
                               upper.post = NA)
  rownames(PaperGamma.tbl) <- c("Equal weights", "maxEnt", "min KL div.", "Hier. prior")
  
  AlphasGamma.tbl <- data.frame(matrix(NA, nrow = 3, ncol = length(av)))
  rownames(AlphasGamma.tbl) <- c("maxEnt", "min KL div.", "Hier. prior")
  colnames(AlphasGamma.tbl) <- paste("alpha", 0:(length(av)-1), sep = "")
    
  ######################################################
  ###### Equal weights alphas
  ######################################################
  alphaEqual <- rep(1/K, K)
  
  ab.Equal.star <- pool.par(alphaEqual, av, bv)
  # Prior
  (PaperGamma.tbl[1, 1:3] <- stat.gamma(ab.Equal.star))
  # Posterior
  (PaperGamma.tbl[1, 4:6] <- stat.gamma(ab.Equal.star + c(ys, n)))
  
  ######################################################
  ##### Maximizing Entropy
  # Maximize H(\pi; alpha),  the entropy of the pooled prior
  ######################################################
  
  optentgamma <- function(alpha, ap, bp){
    entropy.gamma(a = sum(alpha*ap), b = sum(alpha*bp))
  }
  
  optentgamma.inv <- function(alpha.inv, ap, bp){
    alpha <- alpha.01(alpha.inv)
    -optentgamma(alpha, ap, bp)
  }
  
  a.ent <- optim(c(0, 0, 0, 0), optentgamma.inv, ap = av, bp = bv) 
  #            method = "SANN", control=list(maxit = 100000))
  
  alphaMaxEnt.opt <- alpha.01(a.ent$par)
  (AlphasGamma.tbl[1,] <- alphaMaxEnt.opt)
  ab.MaxEnt.star <- pool.par(alphaMaxEnt.opt, av, bv)
  
  # Prior
  (PaperGamma.tbl[2, 1:3] <- stat.gamma(ab.MaxEnt.star))
  
  # Posterior
  (PaperGamma.tbl[2, 4:6] <- stat.gamma(ab.MaxEnt.star + c(ys, n)))
  
  ######################################################
  # Minimizing KL divergence between each density and the pooled prior
  # F = { f_0, f_1, ..., f_K}
  # d_i = KL(f_i || \pi) ; L(alpha) = sum(d_i)
  # minimize the loss function L(F; alpha)
  ######################################################
  
  optklgamma <- function(alpha, ap, bp){
    K <- length(alpha)
    astar <- sum(alpha*ap)
    bstar <- sum(alpha*bp)
    ds <- rep(NA, K) # the distances from each f_i to \pi
    for (i in 1:K){ ds[i] <-  kl.gamma(a0 = ap[i], b0 = bp[i], a1 = astar , b1 = bstar)} 
    return(ds)
  }
  optklgamma.inv <- function(alpha.inv, ap, bp){
    alpha <- alpha.01(alpha.inv)
    sum(optklgamma(alpha, ap, bp))
  }
  
  a.kl <- optim(c(0, 0, 0, 0), optklgamma.inv, ap = av, bp = bv)
  #            method = "SANN", control=list(maxit = 100000))
  
  alphaKL.opt <- alpha.01(a.kl$par)
 (AlphasGamma.tbl[2,] <- alphaKL.opt)
  ab.KL.star <- pool.par(alphaKL.opt, av, bv)
  
  # Prior
  (PaperGamma.tbl[3, 1:3] <- stat.gamma(ab.KL.star))
  
  # Posterior
  (PaperGamma.tbl[3, 4:6] <- stat.gamma(ab.KL.star + c(ys, n)))
  
  ######################################################
  ###### Hierarchical prior
  # \pi(theta|alpha)
  # alpha ~ Dirichlet (X)
  # X = {x_0, x_1, ..., x_K}
  ######################################################
  require("LearnBayes")
  
  M <- 100000
#   X <- c(1, 1, 1, 1, 1)/2 # Jeffrey's prior
  X <- c(1, 1, 1, 1, 1)/length(av)
  alpha.MC <- rdirichlet(M, X)
  gamma.par <- alpha.MC %*% cbind(av, bv)
  
  lambda.par <- apply(gamma.par, 1, function(x) rgamma(1, x[1], x[2]))
  # Prior
  (PaperGamma.tbl[4, 1] <- mean(lambda.par))
  (PaperGamma.tbl[4, 2:3] <- quantile(lambda.par, c(.025, .975)) )
  
  # Posterior 
  library(rstan)
  gammadata.stan <-  list(Y = y, X = X, N = n, K = K, a = av, b = bv)
  # hierpost <- stan(file = "posterior_gamma_pooled.stan",
  #                      data = gammadata.stan, iter = 1, thin = 1, chains = 1)
  # save(hierpost, file = "compiled_gamma_post_sampler.RData")
  load("compiled_gamma_post_sampler.RData")
  hierpostsamp <- stan(fit = hierpost,
                     data = gammadata.stan, iter = 50000, thin = 1, chains = 1)
  posteriors <- extract(hierpostsamp)
  alpha.post <- apply(posteriors$alpha, 2, mean)
  ab.Hier.star <- pool.par(alpha = alpha.post, a = av, b = bv)
  
  (PaperGamma.tbl[4, 4] <- mean(posteriors$lambda))
  (PaperGamma.tbl[4, 5:6] <- quantile(posteriors$lambda, c(.025, .975)) )
  (AlphasGamma.tbl[3,] <- apply(posteriors$alpha, 2, mean) )
  
 ######################################################
  ############## Results
  ######################################################
  # Table
  round(Prior.exps, 3)
  round(PaperGamma.tbl, 3)
  round(AlphasGamma.tbl, 3)
  ###  Plotting
  svg("../plots/gamma_example.svg")
 par(mfrow = c(2, 1))
  # Priors
  curve(fgamma(x, par = ab.Equal.star), 0, 2*lamb, ylab = "Density", main = "Pooled Priors",
        xlab = expression(lambda), lwd = 2)
  curve(fgamma(x, par = ab.MaxEnt.star), 0, 2*lamb, col = 2, add = TRUE, lwd = 2)
  curve(fgamma(x, par = ab.KL.star), 0, 2*lamb, col = 3, add = TRUE, lwd = 2)
  lines(density(lambda.par), col = 4, lwd = 2)
  legend(x = "topright", bty = "n", col = 1:4,
         legend = c("Equal weights (1/K)", "MaxEnt", "MinKL", "Hierarchical"),
         lwd = 2, lty = 1, cex = .8)    
  # Posteriors
  curve(fgamma(x, par = ab.Equal.star + c(ys, n)), 0, 2*lamb, ylab = "Density", main = "Posteriors",
        xlab = expression(lambda), lwd = 2)
  curve(fgamma(x, par = ab.MaxEnt.star + c(ys, n)), 0, 2*lamb, col = 2, add = TRUE, lwd = 2)
  curve(fgamma(x, par = ab.KL.star + c(ys, n)), 0, 2*lamb, col = 3, add = TRUE, lwd = 2)
  lines(density(posteriors$lambda), col = 4, lwd = 2)
  legend(x = "topleft", bty = "n", col = 1:4,
         legend = c("Equal weights (1/K)", "MaxEnt", "MinKL", "Hierarchical"),
         lwd = 2, lty = 1, cex = .8)
  abline( v = ys/n, lwd = 2, lty = 2)
 dev.off()
   
 pars <- list(equal = ab.Equal.star,
              entropy = ab.MaxEnt.star,
              kl = ab.KL.star,
              hier = ab.Hier.star
              )
 
comb.marglikes <- lapply(pars, function(p) ml.gamma(y = y, a = p[1], b = p[2]))
 lapply(comb.marglikes, log)
