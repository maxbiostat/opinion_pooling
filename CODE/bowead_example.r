# P_{t+1} = P_t - C_t + 1.5 * (MSYR) * P_t(1-(P_t/P_0)^2) 
# Inputs: P_0 and MSYR
# Output: P_{1993}
source("pooling_aux.r")
## Prior for P_0 is a Gamma(2.8085, 0.0002886) + 6400, i.e, with offset 6400.
M <- 2E4
P0.prior <- rgamma(M, shape = 2.8085, rate = 0.0002886) + 6400
c(mean(P0.prior), quantile(P0.prior, probs = c(.025, .5, .975)))
summary(P0.prior)
hist(P0.prior, probability = TRUE)

## Prior for MSYR is Gamma(8.2, 372.7)
MSYR.prior <- rgamma(M, shape = 8.2, rate = 372.7)
c(mean(MSYR.prior), quantile(MSYR.prior, probs = c(.025, .5, .975)))
summary(MSYR.prior)
hist(MSYR.prior, probability = TRUE)

## Prior on P_{1993} is N(7800, 1300^2)
P1993.prior <- rnorm(M, mean = 7800, sd = 1300) ## q_2(\phi)
c(mean(P1993.prior), quantile(P1993.prior, probs = c(.025, .5, .975)))
summary(P1993.prior)
hist(P1993.prior, probability = TRUE)
P1993.prior.dens <- density(P1993.prior)

P1993.likelihood2 <- function(x) dnorm(x, m = 8293,  sd = 626) ## L_2(\phi)

CtData <- read.csv("../data/bowhead_kills_1848-2002.csv")
CtData <- subset(CtData, year <= 1992)
barplot(CtData$kills)
getP <- function(msyr, p0, C, Pinit){
  TT <- length(C)
  P <- rep(NA, TT + 1)
  P[1] <- Pinit
  for(t in 1: TT)  P[t+1] <- P[t] - C[t] + 1.5 * msyr * P[t]*(1-(P[t]/p0)^2) 
  return(P)
}
plot(
  seq(min(CtData$year), max(CtData$year) + 1, by = 1),
  getP(msyr = median(MSYR.prior),
       p0 = median(P0.prior),
       C = CtData$kills,
       Pinit = mean(P0.prior)),
  xlab = "year", ylab = "Population size"
)
induced.P1993.prior1 <- unlist(
  parallel::mclapply(1:M, function(i){
    simuP <- getP(msyr = MSYR.prior[i],
                  p0 = P0.prior[i],
                  C = CtData$kills,
                  Pinit = P0.prior[i])
    return(simuP[length(simuP)])
  }, mc.cores = 4)
)
induced.P1993.prior1 <- induced.P1993.prior1[induced.P1993.prior1>0]
induced.P1993.prior1.dens <- density(induced.P1993.prior1)

summary(induced.P1993.prior1)
c(mean(induced.P1993.prior1), quantile(induced.P1993.prior1, probs = c(.025, .5, .975)))
hist(induced.P1993.prior1, probability = TRUE)
distFit <- fGarch::snormFit(induced.P1993.prior1)
curve(fGarch::dsnorm(x, mean = distFit$par["mean"], 
                     sd = distFit$par["sd"], xi =   distFit$par["xi"]), 
      min(induced.P1993.prior1), max(induced.P1993.prior1), col = 2, lwd = 3, add = TRUE)
###########
# Likelihood for ROI
# P1993  = P1978* (1 + ROI)^15
# Likelihood(ROI) = y = g(t8) = exp(a + bt8) -1 where
# t8 is a student t-distribution with df = 8
# g_inverse(y) = t8 = (log(y + 1) - a)/b
# dg_inverse(y)/dy = 1/b(y + 1), hence |J| = abs(1/b(y + 1))
# here a = .0302 and b = .0069
ROI.likelihood <- function(y, a = .0302, b = .0069){
  dt((log(y + 1) - a)/b, df = 8) *  abs(1/(b * (y + 1))) ## L_2(\phi)
} 
ROI.likSamples <- exp(.0302  + .0069 * rt(M, df = 8)) - 1
summary(ROI.likSamples)
c(mean(ROI.likSamples), quantile(ROI.likSamples, probs = c(.025, .5, .975)))

hist(ROI.likSamples, probability = TRUE)
curve(ROI.likelihood, min(ROI.likSamples), max(ROI.likSamples),
      lwd = 3, col = 3, add = TRUE)

# Pzero <- 4765 ## Bowhead pop in 1978 according to table 3 in Brandon & Wade (2008)
Pzero <- 4820 ## Alternatively, we could use table 2 in Givens & Poole (1995) 
induced.P1993.prior2 <- Pzero * (1 + ROI.likSamples)^15
summary(induced.P1993.prior2)
c(mean(induced.P1993.prior2), quantile(induced.P1993.prior2, probs = c(.025, .5, .975)))
hist(induced.P1993.prior2)

induced.P1993.prior2.dens <- density(induced.P1993.prior2)
#########
# The pooled prior LP({q_2(\phi), q_1^\ast(\phi)}, alpha), 
P1993.pooledPrior <- function(x, weights){
  dpoolnorm(x = x,
            D = list(
              f0 = function(x) dnorm(x, mean = 7800, sd = 1300),
              f1 = function(x) fGarch::dsnorm(x, mean = distFit$par["mean"], 
                                              sd = distFit$par["sd"],
                                              xi =   distFit$par["xi"])),
            alphas = weights)
} 
##########
Alpha <- .5
alphas <- c(Alpha, 1-Alpha)
###
plot(P1993.prior.dens, lwd = 3, xlim = c(500, 35000), 
     main = "Pooled distributions", xlab = "Bowhead population in 1993")
curve(fGarch::dsnorm(x, mean = distFit$par["mean"], 
                     sd = distFit$par["sd"], xi = distFit$par["xi"]), 
      500, 35000, col = 2, lwd = 3, add = TRUE)
curve(P1993.pooledPrior(x, weights = c(Alpha, 1-Alpha) ),
      500, 35000,
      col = "blue", lwd = 4, add = TRUE)
legend(x = "topright",
       legend = c("Induced prior q_1*", "prior q_2", "Pooled prior alpha = 1/2"),
                  col = c("black", "red", "blue"), lwd = 2, bty = "n")
#########
# Now let's get the posterior using LP-SpIR [for fixed alpha]
SpIR.posterior <- LP.SpIR(
  k = M,
  l = round(.2 * M),
  Model = function(theta){
    P <- getP(msyr = theta[2], p0 = theta[1],
         C = CtData$kills, Pinit = theta[1])
     return(P[length(P)])
  }, 
  rq1 = function(n) cbind(rgamma(n, shape = 2.8085, rate = 0.0002886) + 6400,
                          rgamma(n, shape = 8.2, rate = 372.7)),
  dq2 = function(x) dnorm(x, m = 7800, sd = 1300),
  dL2 = P1993.likelihood2,
  dL1 = function(x) 1,
  alpha = Alpha
)

# Posterior for the input
Theta.posterior <- SpIR.posterior$theta.resamp
apply(Theta.posterior, 2, function(x) c(mean(x), quantile(x, probs = c(.025, .5, .975))))

c(mean(P0.prior), quantile(P0.prior, probs = c(.025, .5, .975)))
c(mean(MSYR.prior), quantile(MSYR.prior, probs = c(.025, .5, .975)))

# Now the output
P1993.posterior <- SpIR.posterior$phi.resamp
c(mean(P1993.posterior), quantile(P1993.posterior, probs = c(.025, .5, .975)))
c(mean(P1993.prior), quantile(P1993.prior, probs = c(.025, .5, .975)))

hist(P1993.posterior)