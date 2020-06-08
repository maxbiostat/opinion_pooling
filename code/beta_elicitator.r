## Here is a simple method-of-moments based elicitation procedure for the Beta distribution.
## This function will take either 
## (i) prior mean and variance m_0 and v_0, respectively, or;
## (ii) prior mean and coefficient of variation m_0 and c, respectively.
elicit_beta_mean_cv <- function(m0, v0 = NULL, cv = 1){
  if(!is.null(v0)){
    a <- -(m0*v0 + m0^3 - m0^2)/v0
    b <- ((m0-1)*v0 + m0^3 - 2*m0^2 + m0)/v0
  }else{
    a <- -(m0*(cv*m0)^2 + m0^3 - m0^2)/(cv*m0)^2
    b <- ((m0-1)*(cv*m0)^2 + m0^3 - 2*m0^2 + m0)/(cv*m0)^2
  }
  if(a < 0 || b <0) warning("Warning: at least one of the obtained parameters is not valid")
  return(list(a = a, b = b))
}
##
# Testing
# m <- .15
# v <- 0.019
# c0 <- 2
# elicit_beta_mean_cv(m0 = m, v0 = v)
# elicit_beta_mean_cv(m0 = m, cv = c0)

get_parameter_vectors_mean_cv <- function(cv_correct, cv_wrong = .1){
  pars_0 <- elicit_beta_mean_cv(m0 = .1, cv = cv_wrong)
  pars_1 <- elicit_beta_mean_cv(m0 = .2, cv = cv_wrong)
  pars_2 <- elicit_beta_mean_cv(m0 = .5, cv = cv_correct) # big shot
  pars_3 <- elicit_beta_mean_cv(m0 = .8, cv = cv_wrong)
  pars_4 <- elicit_beta_mean_cv(m0 = .9, cv = cv_wrong)
  
  av <- unlist( c(pars_0[1], pars_1[1], pars_2[1], pars_3[1], pars_4[1]) )
  bv <- unlist( c(pars_0[2], pars_1[2], pars_2[2], pars_3[2], pars_4[2]) )
  return(list(
    av = as.numeric(av),
    bv = as.numeric(bv)
  ))
}
#
elicit_beta_median_iq <- function(m, d, q = .90){
  u <- m + d
  loss <- function(x){
    a <- x[1]
    b <- x[2]
    m.hat <- qbeta(.5, shape1 = a, shape2 = b)
    u.hat <- qbeta(q, shape1 = a, shape2 = b)
    error <- .5 * (m.hat - m)^2 + .5*(u.hat-u)^2
    return(error)  
  }
  opt <- suppressWarnings( optim(loss, par = c(1, 1), lower = c(1E-3, 1E-3)) )
  a <- opt$par[1]
  b <- opt$par[2]
  if(a < 0 || b <0) warning("Warning: at least one of the obtained parameters is not valid")
  return(list(a = a, b = b))
}
# mm <- .7
# dd <- .25
# pars <- elicit_beta_median_iq(m = mm, d = dd)
# qbeta(.5, shape1 = pars$a, shape2 = pars$b)
# qbeta(.9, shape1 = pars$a, shape2 = pars$b)
# curve(dbeta(x,shape1 = pars$a, shape2 = pars$b), lwd = 2)
get_parameter_vectors_median_iq <- function(med = 0.5,
                                            iq_correct,
                                            iq_wrong = c(.0125, .025, .05, .1) ){
  pars_0 <- elicit_beta_median_iq(m = med, d = iq_wrong[1])
  pars_1 <- elicit_beta_median_iq(m = med, d = iq_wrong[2])
  pars_2 <- elicit_beta_median_iq(m = med, d = iq_correct) # big shot
  pars_3 <- elicit_beta_median_iq(m = med, d = iq_wrong[3])
  pars_4 <- elicit_beta_median_iq(m = med, d = iq_wrong[4])
  
  av <- unlist( c(pars_0[1], pars_1[1], pars_2[1], pars_3[1], pars_4[1]) )
  bv <- unlist( c(pars_0[2], pars_1[2], pars_2[2], pars_3[2], pars_4[2]) )
  return(list(
    av = as.numeric(av),
    bv = as.numeric(bv)
  ))
}
#
plot_densities <- function(pars, lg = TRUE){
  av <- pars$av
  bv <- pars$bv
  K <- length(av)
  curve(dbeta(x, av[1], bv[1]), lwd = 2)
  for(k in 2:K){
    curve(dbeta(x, av[k], bv[k]), lwd = 2, col = k, add = TRUE)
  }
  if(lg){
    legend(x = "top", col = 1:K, lwd = 2, lty = 1, bty = 'n', legend = paste("Expert_", (1:K)-1, sep = "") )
  }
}