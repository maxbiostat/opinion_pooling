# Leo Bastos & Luiz Max Carvalho (2017)
source("pooling_aux.r")
source("beta_elicitator.r")
get_odds <- function(x) {
  y <- x[order(x, decreasing = TRUE)]
  y[1]/y[2]
}
beta.sd <- function(a, b) sqrt( (a*b)/( (a + b + 1) * (a + b)^2 ))

compute_odds <- function(av, bv, X, y, n){
  K <- length(av)
 require(rstan)
  betadata.stan <-  list(Y = y, X = X, N = n, K = K, a = av, b = bv)
  betadata.stan.exp <-  list(Y = y,
                             means = digamma(X)-digamma(X[K]),
                             Sigma = constructSigma(X),
                             N = n, K = K, a = av, b = bv)
  
  load("compiled_beta_post_Dirichlet_sampler.RData")
  load("compiled_beta_post_logisticNormal_sampler.RData")
  capture.output( hierpostsamp.dir <- stan(fit= hierpost.dir,
                           data = betadata.stan, iter = 5000, thin = 1, chains = 1) 
                  )
  capture.output( hierpostsamp.exp <- stan(fit= hierpost.exp,
                           data = betadata.stan.exp, iter = 5000, thin = 1, chains = 1)
                  )
  
  posteriors.dir <- extract(hierpostsamp.dir)
  posteriors.exp <- extract(hierpostsamp.exp)
  
  alphas.exp <- matrix(NA, nrow = nrow(posteriors.exp$m), ncol = K )
  for (i in 1: nrow(posteriors.exp$m)){
    alphas.exp[i, ] <- exp(posteriors.exp$m[i, ])/ sum(exp(posteriors.exp$m[i, ])) 
  }
  # post.alpha.dir.cred <- apply(posteriors.dir$alpha, 2, quantile, probs = c(.025, .975))
  # post.alpha.exp.cred <- apply(alphas.exp, 2, quantile, probs = c(.025, .975))
  
  post.alpha.dir <- apply(posteriors.dir$alpha, 2, mean)
  post.alpha.exp <- apply(alphas.exp, 2, mean)
  return(
    data.frame(
      dirichlet = get_odds(post.alpha.dir),
      aitchinson = get_odds(post.alpha.exp)
    )
  )
}

########################################
########################################
p0 <- elicit.beta(m0 = 4/10, v0 = 1/1000); a0 <- p0$a ; b0 <- p0$b
p1 <- elicit.beta(m0 = 4/10, v0 = 1/100); a1 <- p1$a ; b1 <- p1$b
p2 <- elicit.beta(m0 = 4/10, v0 = 1/50); a2 <- p2$a ; b2 <- p2$b
p3 <- elicit.beta(m0 = 4/10, v0 = 1/10); a3 <- p3$a ; b3 <- p3$b

As <- c(a0, a1, a2, a3)
Bs <- c(b0, b1, b2, b3)
( SDs <- beta.sd(As, Bs) )

DirPar <- c(1, 1, 1, 1)/2 # Jeffreys' prior

# Data
count <- 4
total <- 10

Nrep <- 100
system.time(
  suppressWarnings(
  SimusList <-mclapply(1:Nrep,
                       function(i){
                         set.seed(i)
                         compute_odds(av = As, bv = Bs, X = DirPar, y = count, n = total)
                       },
                       mc.cores = 10)
  )
)

SimuDt <- data.table::rbindlist(SimusList)

plot(SimuDt)
apply(SimuDt, 2, mean)
apply(SimuDt, 2, sd)
apply(SimuDt, 2, quantile, probs = c(.025, .975))

