logit <- function(p) log(p)-log(1-p)
inv_logit <- function(l) exp(l)/(1 + exp(l))
dbeta_logit <- function(x, a, b, log = FALSE){
  ans <- dbeta(inv_logit(x), a, b, log = TRUE) + x - 2*log1p(exp(x))
  if(!log) ans <- exp(ans)
  return(ans)
}
dbeta_logit <- Vectorize(dbeta_logit)
#############
a <- 8
b <- 1
X <- rbeta(1e6, a , b)
Y <- logit(X)

hist(Y, probability = TRUE)
curve(dbeta_logit(x, a = a, b = b), min(Y), max(Y), lwd = 2, add = TRUE)
