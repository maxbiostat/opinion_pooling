source("pooling_aux.r")

### Functions
getP <- function(msyr, p0, C, Pinit){
  TT <- length(C)
  P <- rep(NA, TT + 1)
  P[1] <- Pinit
  for(t in 1: TT)  P[t+1] <- P[t] - C[t] + 1.5 * msyr * P[t]*(1-(P[t]/p0)^2) 
  return(P)
}
#
get_P1993 <- function(theta){
  P <- getP(msyr = theta[2], p0 = theta[1],
            C = CtData$kills, Pinit = theta[1])
  return(P[length(P)])
}
#
P1993_likelihood2 <- function(x) dnorm(x, mean = 8293,  sd = 626) ## L_2(\phi)
#
ROI_likelihood <- function(y, a = .0302, b = .0069, log = FALSE){
  if(is.nan(y) && log == FALSE) return(0)
  if(is.nan(y) && log == TRUE) return(-Inf)
  ## L_2(\phi)
  if(log){
    dt((log(y + 1) - a)/b, df = 8, log = TRUE) + log(abs(1/(b * (y + 1)))) 
  }else{
    dt((log(y + 1) - a)/b, df = 8) *  abs(1/(b * (y + 1))) 
  }
}
#
get_ROI <- function(pars){
  Ps <-  getP(msyr = pars[2], p0 = pars[1],
              C = CtData$kills, Pinit = pars[1])
  P1978 <- Ps[length(Ps) - 15]
  P1993 <- Ps[length(Ps)]
  ROI <- (P1993/P1978)^{1/15} - 1
  return(ROI)
}
get_ROI <- compiler::cmpfun(get_ROI)
#################### Analysis
CtData <- read.csv("../data/bowhead_kills_1848-2002.csv")
CtData <- subset(CtData, year <= 1992)
barplot(CtData$kills)
M <- 1E5
SpIR.posterior <- LPSpIR_varying_alpha(
  k = M,
  l = M,
  Model = get_P1993, 
  rq1 = function(n) cbind(rgamma(n, shape = 2.8085, rate = 0.0002886) + 6400,
                          rgamma(n, shape = 8.2, rate = 372.7),
                          rbeta(n, shape1 = 1, shape2 = 1)),
  dq2 = function(x) dnorm(x, m = 7800, sd = 1300),
  dL2 = P1993_likelihood2,
  dL1 = function(x) ROI_likelihood(get_ROI(x[1:2])),
  cores = 6
)

Theta.posterior <- data.frame(SpIR.posterior$theta.resamp)
names(Theta.posterior) <- c("P0", "MSYR", "alpha")
Theta.posterior <- data.frame(
  Theta.posterior,
  ROI = apply(Theta.posterior, 1, function(x) get_ROI(x[1:2]))
)

tempAlpha <- Theta.posterior$alpha
Theta.posterior$alpha <- NULL
Theta.posterior$P1993 <- SpIR.posterior$phi.resamp
Theta.posterior$alpha <- tempAlpha

MCMC.theta <- coda::as.mcmc(Theta.posterior)
plot(MCMC.theta)
summary(MCMC.theta)

Theta.posterior$alpha_status <- "varying"
Theta.posterior$method <- "SpIR"

write.csv(Theta.posterior, 
          file = "../data/output/bowhead_results_SpIR_varying.csv",
          row.names = FALSE)
