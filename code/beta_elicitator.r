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
######
elicit_beta_median_iq <- function(m, d, q = .90){
  u <- m + d
  loss <- function(x){
    a <- x[1]
    b <- x[2]
    m.hat <- qbeta(.5, shape1 = a, shape2 = b)
    u.hat <- qbeta(q, shape1 = a, shape2 = b)
    error <-  (m.hat - m)^2 + (u.hat-u)^2
    return(error)  
  }
  opt <- suppressWarnings( optim(loss, par = c(1, 1), lower = c(1E-3, 1E-3)) )
  a <- opt$par[1]
  b <- opt$par[2]
  if(a < 0 || b < 0) warning("Warning: at least one of the obtained parameters is not valid")
  return(list(a = a, b = b))
}
### TEST
# mm <- .7
# dd <- .25
# pars <- elicit_beta_median_iq(m = mm, d = dd)
# qbeta(.5, shape1 = pars$a, shape2 = pars$b)
# qbeta(.9, shape1 = pars$a, shape2 = pars$b)
# curve(dbeta(x,shape1 = pars$a, shape2 = pars$b), lwd = 2)

######
elicit_beta_quantiles <- function(l, u, alpha = .95){
  q0 <- (1-alpha)/2
  q1 <- (1 + alpha)/2
  loss <- function(x){
    a <- x[1]
    b <- x[2]
    l.hat <- qbeta(q0, shape1 = a, shape2 = b)
    u.hat <- qbeta(q1, shape1 = a, shape2 = b)
    # error <- (l.hat - l)^2 + (u.hat-u)^2
    error <- abs(l.hat - l) + abs(u.hat-u) ## L1 norm: better?
    return(error)  
  }
  opt <- suppressWarnings( optim(loss, par = c(1, 1), lower = c(1E-3, 1E-3)) )
  a <- opt$par[1]
  b <- opt$par[2]
  if(a < 0 || b < 0) warning("Warning: at least one of the obtained parameters is not valid")
  return(list(a = a, b = b))
}
## TEST
# alpha <- .99
# ll <- .15
# uu <- .2
# pars <- elicit_beta_quantiles(l = ll, u = uu, alpha = alpha)
# pars$a / (pars$a + pars$b) # mean
# qbeta((1-alpha)/2, shape1 = pars$a, shape2 = pars$b)
# qbeta(.5, shape1 = pars$a, shape2 = pars$b)
# qbeta((1 + alpha)/2, shape1 = pars$a, shape2 = pars$b)
# curve(dbeta(x,shape1 = pars$a, shape2 = pars$b), lwd = 2)

