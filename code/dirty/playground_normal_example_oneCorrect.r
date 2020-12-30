##################
source("../pooling_aux.r")

set.seed(666)
n <- 10
truemean <- 3
fixed.sigma <- 1
y.obs <- rnorm(n = n, mean = truemean, sd = fixed.sigma)

m0s <- c(1, 2, 3, 4, 5)


r_cv <- function(c2, log = FALSE){
  cvs <- rep(.1, length(m0s))
  cvs[3] <- c2
  v0s <-(cvs*m0s)^2
  golden.weights <- get_normal_weights(obs = y.obs,
                                       sigma = fixed.sigma, ms = m0s, vs = v0s)
  x <- golden.weights$mls
  y <- x[-3][order(x[-3], decreasing = TRUE)]
  r <- x[3]-y[1]
  if(!log){
   r <- exp(r)
  }
  return(r)
}
r_cv <- Vectorize(r_cv)

curve(r_cv, 0, 3, 
      ylab = "correct_expert/second_best",
      xlab = "Coefficient of variation")


# l_cv <- function(cv, log = FALSE){
#    v <- (cv*m0s[3])^2
#    ans <- normal_mean_marg_like(y = y.obs, sigma = fixed.sigma, m = m0s[3],
#                                v = v, log = TRUE)
#    if(!log) ans <- exp(ans)
#    return(ans)
# }
# l_cv <- Vectorize(l_cv)
# 
# curve(l_cv(x, log = TRUE), 0, 3)
# pool_par_gauss(alpha = golden.weights$weights, m = m0s, v = v0s)
