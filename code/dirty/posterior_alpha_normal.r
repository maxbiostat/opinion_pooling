compute_I <- function(y, sigmasq, m0, v0, log = FALSE){
  xbar <- mean(y)
  s2 <-  sum((y-xbar)^2)
  vtilde <- n*v0 + sigmasq 
  l1 <- -s2/(2*sigmasq)
  l2 <- - (n*(xbar-m0)^2)/(2*vtilde)
  l3 <- -0.5 * ( log(vtilde) + (n-1) * log(2*pi*sigmasq) + log(2*pi))
  ans <- l1 + l2 + l3
  if(!log) ans <- exp(ans)
  return(ans)
}
##
truemu <-  2
sigmasq <- 1
n <- 20
set.seed(666)
y <- rnorm(n = n, mean = truemu, sd = sqrt(sigmasq))

m0 <- 0
v0 <- 1

post_kernel <- function(mu){
  ldens <- dnorm(mu, mean = m0, sd = sqrt(v0), log = TRUE) + 
    sum(dnorm(y, mean = mu, sd = sqrt(sigmasq), log = TRUE))
    return(exp(ldens))
}
post_kernel <- Vectorize(post_kernel)

curve(post_kernel, truemu-3, truemu + 3)

integrate(post_kernel, -Inf, Inf, subdivisions = 1E5L)
compute_I(y = y, sigmasq = sigmasq, m0 = m0, v0 = v0)
compute_I(y = y, sigmasq = sigmasq, m0 = m0, v0 = v0, log = TRUE)
