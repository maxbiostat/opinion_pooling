## Here is a simple method-of-moments based elicitation procedure for the Beta distribution.
## This function will take either 
## (i) prior mean and variance m_0 and v_0, respectively, or;
## (ii) prior mean and coefficient of variation m_0 and c, respectively.
elicit_beta <- function(m0, v0 = NULL, cv = 1){
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
# elicit_beta(m0 = m, v0 = v)
# elicit_beta(m0 = m, cv = c0)

get_parameter_vectors <- function(cv_correct, cv_wrong = .1){
  pars_0 <- elicit_beta(m0 = .1, cv = cv_wrong)
  pars_1 <- elicit_beta(m0 = .2, cv = cv_wrong)
  pars_2 <- elicit_beta(m0 = .5, cv = cv_correct) # big shot
  pars_3 <- elicit_beta(m0 = .8, cv = cv_wrong)
  pars_4 <- elicit_beta(m0 = .9, cv = cv_wrong)
  
  av <- unlist( c(pars_0[1], pars_1[1], pars_2[1], pars_3[1], pars_4[1]) )
  bv <- unlist( c(pars_0[2], pars_1[2], pars_2[2], pars_3[2], pars_4[2]) )
  return(list(
    av = as.numeric(av),
    bv = as.numeric(bv)
  ))
}

plot_densities <- function(pars){
  av <- pars$av
  bv <- pars$bv
  K <- length(av)
  curve(dbeta(x, av[1], bv[1]), lwd = 2)
  for(k in 2:K){
    curve(dbeta(x, av[k], bv[k]), lwd = 2, col = k, add = TRUE)
  }
}