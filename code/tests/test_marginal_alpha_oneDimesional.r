source("pooling_aux.r")
source("beta_elicitator.r")

p0 <- elicit_beta_mean_cv(m0 = 6/10, cv = .5); a0 <- p0$a ; b0 <- p0$b
p1 <- elicit_beta_mean_cv(m0 = 6/10, cv = .25); a1 <- p1$a ; b1 <- p1$b

av <- c(a0, a1)
bv <- c(b0, b1)
K <- length(av)

n <- 10
y <- 9

core <- function(a){
  alpha <- c(a, 1-a)
  ppars <- pool_par(alpha, a = av, b = bv)
  ans <- lbeta(a = ppars[1] + y, b = ppars[2] + n-y) - lbeta(a = ppars[1], b = ppars[2]) + lchoose(n, y)
  ans <- exp(ans)
  return(ans)
}
core <- Vectorize(core)
curve(core)

M <- 1e5
( Z <- integrate(core, 0, 1)$value )
A.samples <- rbeta(M, 1, 1)
mean(core(A.samples))


integrate(function(x) x * core(x)/Z, 0, 1)

AA.samples <- cbind(A.samples, 1-A.samples)

core2 <- function(alpha){
  ppars <- pool_par(alpha, a = av, b = bv)
  ans <- lbeta(a = ppars[1] + y, b = ppars[2] + n-y) - lbeta(a = ppars[1], b = ppars[2]) + lchoose(n, y)
  ans <- exp(ans)
  return(ans)
}

densities <- apply(AA.samples, 1, core2)

( lZ <- log(mean(densities)) )
log(Z)


phi <- function(alpha){
  ppars <- pool_par(alpha, a = av, b = bv)
  ans <- log(alpha) + lbeta(a = ppars[1] + y, b = ppars[2] + n-y) - lbeta(a = ppars[1], b = ppars[2]) + lchoose(n, y) - lZ
  ans <- exp(ans)
  return(ans)
}

rowMeans(apply(AA.samples, 1, phi))


