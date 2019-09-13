# P_{t+1} = P_t - C_t + 1.5 * (MSYR) * P_t(1-(P_t/P_0)^2) 
# Inputs: P_0 and MSYR
# Output: P_{1993}
# Obtain the posterior via Sampling-importance-resampling (SpIR)
source("pooling_aux.r")
CtData <- read.csv("../data/bowhead_kills_1848-2002.csv")
CtData <- subset(CtData, year <= 1992)

getPosterior <- function(M, Alpha, cores = 4){
  require(compiler)
  cat("Doing alpha = ", Alpha, "\n")
  Result <- vector(0, mode = "list")
  
  P1993_likelihood2 <- function(x) dnorm(x, m = 8293,  sd = 626) ## L_2(\phi)
  P1993_likelihood2 <- cmpfun(P1993_likelihood2)
  #
  get_ROI <- function(pars){
    Ps <-  getP(msyr = pars[2], p0 = pars[1],
                C = CtData$kills, Pinit = pars[1])
    P1978 <- Ps[length(Ps) - 15]
    P1993 <- Ps[length(Ps)]
    ROI <- (P1993/P1978)^{1/15} - 1
    return(ROI)
  }
  get_ROI <- cmpfun(get_ROI)
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
  ROI_likelihood <- cmpfun(ROI_likelihood)
  #
  getP <- function(msyr, p0, C, Pinit){
    TT <- length(C)
    P <- rep(NA, TT + 1)
    P[1] <- Pinit
    for(t in 1: TT)  P[t+1] <- P[t] - C[t] + 1.5 * msyr * P[t]*(1-(P[t]/p0)^2) 
    return(P)
  }
  getP <- cmpfun(getP)
  #
  get_P1993 <- function(theta){
    P <- getP(msyr = theta[2], p0 = theta[1],
              C = CtData$kills, Pinit = theta[1])
    return(P[length(P)])
  }
  get_P1993 <- cmpfun(get_P1993)
  #
  SpIR.posterior <- LP_SpIR(
    k = M,
    l = round(.2 * M),
    Model = get_P1993, 
    rq1 = function(n) cbind(rgamma(n, shape = 2.8085, rate = 0.0002886) + 6400,
                            rgamma(n, shape = 8.2, rate = 372.7)),
    dq2 = function(x) dnorm(x, m = 7800, sd = 1300),
    dL2 = P1993_likelihood2,
    dL1 = function(x) ROI_likelihood(get_ROI(x)),
    alpha = Alpha,
    cores = cores
  )
  # Posterior for the input
  Theta.posterior <- data.frame(SpIR.posterior$theta.resamp)
  names(Theta.posterior) <- c("P0", "MSYR")
  Theta.posterior <- data.frame(
    Theta.posterior,
    ROI = apply(Theta.posterior, 1, get_ROI)
  )
  #
  Posterior <- data.frame(Theta.posterior,
                          P1993 = SpIR.posterior$phi.resamp) 
  Result$samples <- Posterior
  Result$quantiles <- apply(Posterior, 2,
                            function(x) c(mean = mean(x), quantile(x, probs = c(.025, .5, .975))))
  
  Result$correlation <- cor(Posterior)
  return(Result)
}

alphas <- c(0, .25,  .5, .75, 1)
# For some reason the code below throws a really weird error. Can't be
simu.time <- system.time(
  Simus <- lapply(alphas, function(x) getPosterior(M = 1e5, Alpha = x, cores = 6)  )
)

simu.time
#########
Samples <- lapply(Simus, function(x) x$samples)
for(i in seq_along(alphas)) Samples[[i]]$alpha <- alphas[i]

library(reshape2)
forplot <- melt( do.call(rbind, Samples), id.vars = "alpha", variable.name = "parameter")
forplot$alpha <- as.factor(forplot$alpha)

library(ggplot2)

ggplot(forplot, aes (x = value, colour = alpha, fill = alpha)) +
  geom_density(alpha = .4) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  facet_wrap(.~parameter, scales = "free")+
  theme_bw(base_size = 16)

alpha_sens_boxplots <- ggplot(forplot, aes (x = alpha, y = value, colour = alpha, fill = alpha)) +
  geom_boxplot(alpha = .4) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_discrete(expression(bold(alpha)), expand = c(0, 0)) +
  facet_wrap(.~parameter, scales = "free")+
  theme_bw(base_size = 16)
alpha_sens_boxplots

ggsave(plot = alpha_sens_boxplots, file = "../plots/alpha_sensitivity_bowhead.pdf")
