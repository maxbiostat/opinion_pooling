## Here is a simple method-of-moments based elicitation procedure for the Beta distribution.
## This function will take either 
## (i) prior mean and variance m_0 and v_0, respectively, or;
## (ii) prior mean and coefficient of variation m_0 and c, respectively.
elicit.beta <- function(m0, v0 = NULL, c = 1){
  m <- (m0-1)/m0
  if(!is.null(v0)){
    a <- -(m + v0*(1-m)^2)/(v0*(1-m)^3)
    b <- -m*(-(m + v0*(1-m)^2)/(v0*(1-m)^3))
  }else{
    a <- -(m + (c*m0)^2*(1-m)^2)/((c*m0)^2*(1-m)^3)
    b <- -m*(-(m + (c*m0)^2*(1-m)^2)/((c*m0)^2*(1-m)^3))
  }
  return(list(a = a, b = b))
}
##
# Testing
#M <- .15
#V <- .019
#c0 <- 2
#elicit.beta(m0 = M, v0 = V)
#elicit.beta(m0 = M, c = c0)
