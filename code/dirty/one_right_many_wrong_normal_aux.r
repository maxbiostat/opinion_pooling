###################################################
get_parameter_vectors_mean_cv <- function(cv_correct, cv_wrong = .1){
  mv <- c(1, 2, 3, 4, 5)
  cvs <- rep(.1, length(mv))
  cvs[3] <- cv_correct
  v0s <-(cvs*mv)^2
  return(list(
    mv = as.numeric(mv),
    vv = as.numeric(v0s)
  ))
}
#
plot_densities <- function(pars, lg = TRUE){
  mv <- pars$mv
  vv <- pars$vv
  xmin <- max(mv) - 2*max(mv) 
  xmax <- max(mv) + 2*max(mv)  
  K <- length(mv)
  curve(dnorm(x, mean = mv[1],  sd = sqrt(vv[1]) ), xlim = c(xmin, xmax), lwd = 2)
  for(k in 2:K){
    curve(dnorm(x, mean = mv[k], sd = sqrt(vv[k])), lwd = 2, col = k, add = TRUE)
  }
  if(lg){
    legend(x = "top", col = 1:K, lwd = 2, lty = 1, bty = 'n', legend = paste("Expert_", (1:K)-1, sep = "") )
  }
}
#