## Here is a simple method-of-moments based elicitation procedure for the Beta distribution.
## This function will take either 
## (i) prior mean and variance m_0 and v_0, respectively, or;
## (ii) prior mean and coefficient of variation m_0 and c, respectively.
elicit.beta <- function(m0, v0 = NULL, cv = 1){
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
m <- .15
v <- 0.019
c0 <- 2
elicit.beta(m0 = m, v0 = v)
elicit.beta(m0 = m, cv = c0)