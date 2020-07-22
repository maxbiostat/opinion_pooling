f2_unnorm <- function(x, log = FALSE){
  ans <- log(alphas[2])  -(x-2.5)^2/(2 * 2^2) 
  if(!log) ans <- exp(ans)
  return(ans)
} 
f2_unnorm <- Vectorize(f2_unnorm)
m1 <- simple_is(samples = X1, target_dens = f2_unnorm, prop_dens = f1)$mu
m3 <- simple_is(samples = X3, target_dens = f2_unnorm, prop_dens = f3)$mu
m4 <- simple_is(samples = X4, target_dens = f2_unnorm, prop_dens = f4)$mu
m1; m3; m4
integrate(f2_unnorm, 0, Inf)
