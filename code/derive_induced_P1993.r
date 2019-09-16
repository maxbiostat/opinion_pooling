# P_{t+1} = P_t - C_t + 1.5 * (MSYR) * P_t(1-(P_t/P_0)^2) 
# Inputs: P_0 and MSYR
# Output: P_{1993}
source("pooling_aux.r")
## Prior for P_0 is a Gamma(2.8085, 0.0002886) + 6400, i.e, with offset 6400.
M <- 1E5
P0.prior <- rgamma(M, shape = 2.8085, rate = 0.0002886) + 6400
hist(P0.prior, probability = TRUE)

## Prior for MSYR is Gamma(8.2, 372.7)
MSYR.prior <- rgamma(M, shape = 8.2, rate = 372.7)
hist(MSYR.prior, probability = TRUE)

## Prior on P_{1993} is N(7800, 1300^2)
P1993.prior <- rnorm(M, mean = 7800, sd = 1300) ## q_2(\phi)
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


## plot
hist(induced.P1993.prior1, probability = TRUE, xlab = expression(P[1993]), main = "")
distFit <- fGarch::snormFit(induced.P1993.prior1)
curve(fGarch::dsnorm(x, mean = distFit$par["mean"], 
                     sd = distFit$par["sd"], xi =   distFit$par["xi"]), 
      min(induced.P1993.prior1), max(induced.P1993.prior1), col = 2, lwd = 3, add = TRUE)
lnfit <- fitdistrplus::fitdist(induced.P1993.prior1, dlnorm)
curve(dlnorm(x, meanlog = lnfit$estimate["meanlog"], sdlog = lnfit$estimate["sdlog"] ),
      min(induced.P1993.prior1), max(induced.P1993.prior1), col = 3, lwd = 3, add = TRUE)
nfit <- fitdistrplus::fitdist(induced.P1993.prior1, dnorm)
curve(dnorm(x, mean = nfit$estimate["mean"], sd = nfit$estimate["sd"] ),
      min(induced.P1993.prior1), max(induced.P1993.prior1), col = 4, lwd = 3, add = TRUE)
legend(x = "topright",
       legend = c("Normal", "Log-normal", "Skew-normal"), lwd = 2, bty = 'n', col = c(4, 3, 2))

### Looking at the fit at the tails
alpha <- .975
qq <- quantile(induced.P1993.prior1, probs = alpha)
fGarch::psnorm(qq, mean = distFit$par["mean"], sd = distFit$par["sd"], xi =   distFit$par["xi"])
plnorm(qq, meanlog = lnfit$estimate["meanlog"], sdlog = lnfit$estimate["sdlog"] )
pnorm(qq, mean = nfit$estimate["mean"], sd = nfit$estimate["sd"] )

## AIC normal
k_normal <- 2
2 * k_normal - 2*nfit$loglik # could simply do nfit$AIC
## AIC skew-normal
k_snormal <- 3
2 * k_snormal - 2*-distFit$objective
## AIC 
k_lnormal <- k_normal
2 * k_lnormal - 2*lnfit$loglik 
## conclusion is that while the skew normal is probably a better fit, the normal distribution is quite good and allows for more analytical work. 
