CtData <- read.csv("../data/bowhead_kills_1848-2002.csv")
CtData <- subset(CtData, year <= 1992)
stanData <- list(
  p0_shape = 2.8085,
  p0_rate = .0002886,
  msyr_shape = 8.2,
  msyr_rate = 372.7,
  a = .0302,
  b = .0069,
  roi_df = 8,
  phi_mu = 7800,
  phi_sd = 1300,
  ind_mu = 18289.03,
  ind_sd = 6085.33,
  ind_xi = 1.12,
  phi_lik_mu = 8293, 
  phi_lik_sd = 626,
  N = length(CtData$kills),
  C = CtData$kills,
  alpha = 1/2
)
library(rstan)
Nit <- 2E5
# fixed_alpha <- stan(file = "bowhead.stan",
#                     data = stanData, iter = 1,
#                     pars = c("P0", "MSYR", "ROI", "P1993"),
#                     init = list(c1 = list(P0 = 14000, MSYR = .021)),
#                     thin = 1, chains = 1)
# save(fixed_alpha, file = "compiled_bowhead_fixedAlpha.RData")
load("compiled_bowhead_fixedAlpha.RData")
fixed_alpha_sampling <- stan(fit = fixed_alpha, data = stanData,
                             pars = c("P0", "MSYR", "ROI", "P1993"),
                             init = list(c1 = list(P0 = 18000, MSYR = .011),
                                         c2 = list(P0 = 10000, MSYR = .021)),
                             control = list(adapt_delta = .95),
                             iter = Nit, thin = 1, chains = 2)
fixedAlpha_posterior <- as.data.frame(extract(fixed_alpha_sampling))
apply(fixedAlpha_posterior[, 1:4], 2, coda::effectiveSize)
apply(fixedAlpha_posterior[, 1:4], 2,
      function(x) c(mean = mean(x), quantile(x, probs = c(.025, .5, .975))))
coda::traceplot(As.mcmc.list(fixed_alpha_sampling))
plot(fixedAlpha_posterior[, 1:4])
#####################
# Now sampling with alpha
stanData2 <- stanData
stanData2$alpha <- NULL
stanData2$a_alpha <- 1
stanData2$b_alpha <- 1
# bowhead_dirichlet <- stan(file = "bowhead_alpha_dirichlet.stan",
#                     data = stanData2, iter = 1,
#                     pars = c("P0", "MSYR", "ROI", "P1993", "alpha"),
#                     init = list(c1 = list(P0 = 14000, MSYR = .021, alpha = 1/2)),
#                     thin = 1, chains = 1)
# save(bowhead_dirichlet, file = "compiled_bowhead_Dirichilet.RData")
load("compiled_bowhead_Dirichilet.RData")
dirichlet_sampling <- stan(fit = bowhead_dirichlet, data = stanData2,
                           pars = c("P0", "MSYR", "ROI", "P1993", "alpha"),
                             init = list(c1 = list(P0 = 18000, MSYR = .011, alpha = 1/2),
                                         c2 = list(P0 = 10000, MSYR = .021, alpha = 1/2)),
                           control = list(adapt_delta = .85),
                             iter = Nit, thin = 1, chains = 2)
Dirichlet_posterior <- as.data.frame(extract(dirichlet_sampling))
apply(Dirichlet_posterior[, 1:5], 2, coda::effectiveSize)
apply(Dirichlet_posterior[, 1:5], 2,
      function(x) c(mean = mean(x), quantile(x, probs = c(.025, .5, .975))))
coda::traceplot(As.mcmc.list(dirichlet_sampling))
plot(Dirichlet_posterior[, 1:5])