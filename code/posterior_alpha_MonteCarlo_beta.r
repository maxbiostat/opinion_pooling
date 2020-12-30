source("pooling_aux.r")
source("beta_elicitator.r")

p0 <- elicit_beta_mean_cv(m0 = 5/10, cv = .5); a0 <- p0$a ; b0 <- p0$b
p1 <- elicit_beta_mean_cv(m0 = 5/10, cv = .25); a1 <- p1$a ; b1 <- p1$b
p2 <- elicit_beta_mean_cv(m0 = 5/10, cv = .125); a2 <- p2$a ; b2 <- p2$b
p3 <- elicit_beta_mean_cv(m0 = 5/10, cv = .0625); a3 <- p3$a ; b3 <- p3$b

av <- c(a0, a1, a2, a3)
bv <- c(b0, b1, b2, b3)
K <- length(av)

n <- 100
y <- 50

require("LearnBayes")

M <- 10000
X <- c(1, 1, 1, 1)/20

set.seed(666)

alpha.MC.dirichlet <- rdirichlet(M, X)
alpha.MC.logisticNormal <- rlogisticnorm(N = M,
                                         m = digamma(X)-digamma(X[K]),
                                         Sigma = constructSigma(X))

apply(alpha.MC.dirichlet, 2, mean)
apply(alpha.MC.logisticNormal, 2, mean)

apply(alpha.MC.dirichlet, 2, sd)
apply(alpha.MC.logisticNormal, 2, sd)

logkernel <- function(alpha, data){
  y <- data[1]
  n <- data[2]
  ppars <- pool_par(alpha, a = av, b = bv)
  ans <- lbeta(a = ppars[1] + y, b = ppars[2] + n-y) - lbeta(a = ppars[1], b = ppars[2])
  return(ans)
}

lZ.dirichlet <- matrixStats::logSumExp(apply(alpha.MC.dirichlet, 1, logkernel, data = c(y, n)))-log(M) 
lZ.logisticNormal <- matrixStats::logSumExp(apply(alpha.MC.logisticNormal, 1, logkernel, data = c(y, n)))-log(M) 


phi <- function(alpha, data, lZ){
  y <- data[1]
  n <- data[2]
  ppars <- pool_par(alpha, a = av, b = bv)
  ans <- log(alpha) + lbeta(a = ppars[1] + y, b = ppars[2] + n-y) - lbeta(a = ppars[1], b = ppars[2]) - lZ
  ans <- exp(ans)
  return(ans)
}

round(post.phi.dirichlet <- rowMeans(apply(alpha.MC.dirichlet, 1, phi, data = c(y, n), lZ  = lZ.dirichlet)), 3 )
round(post.phi.logisticNormal <- rowMeans(apply(alpha.MC.logisticNormal, 1, phi, data = c(y, n), lZ = lZ.logisticNormal)), 3)

