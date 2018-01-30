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
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Nit <- 2E5
fixed_alpha <- stan(file = "../code/stan/bowhead.stan",
                    data = stanData, iter = 1,
                    pars = c("P0", "MSYR", "ROI", "P1993"),
                    init = list(c1 = list(P0 = 18000, MSYR = .011)), 
                    thin = 1, chains = 1)
# save(fixed_alpha, file = "../code/stan/compiled_bowhead_fixedAlpha.RData")
# load("../code/stan/compiled_bowhead_fixedAlpha.RData")
fixed_alpha_sampling <- stan(fit = fixed_alpha, data = stanData,
                             pars = c("P0", "MSYR", "ROI", "P1993"),
                             init = list(c1 = list(P0 = 17000, MSYR = .021),
                                         c2 = list(P0 = 16000, MSYR = .017),
                                         c3 = list(P0 = 16000, MSYR = .020),
                                         c4 = list(P0 = 18000, MSYR = .015)),
                             control = list(adapt_delta = .95),
                             iter = Nit, thin = 1, chains = 4)
fixed_alpha_sampling
stan_trace(fixed_alpha_sampling)
pairs(fixed_alpha_sampling, pars = c("P0", "MSYR", "ROI", "P1993"))

#####################
# Now sampling with alpha
stanData2 <- stanData
stanData2$alpha <- NULL
stanData2$a_alpha <- 1
stanData2$b_alpha <- 1
bowhead_dirichlet <- stan(file = "../code/stan/bowhead_alpha_dirichlet.stan",
                    data = stanData2, iter = 1,
                    pars = c("P0", "MSYR", "ROI", "P1993", "alpha"),
                    init = list(c1 = list(P0 = 14000, MSYR = .021, alpha = 1/2)),
                    thin = 1, chains = 1)
# save(bowhead_dirichlet, file = "../code/stan/compiled_bowhead_Dirichilet.RData")
# load("../code/stan/compiled_bowhead_Dirichilet.RData")
dirichlet_sampling <- stan(fit = bowhead_dirichlet, data = stanData2,
                           pars = c("P0", "MSYR", "ROI", "P1993", "alpha"),
                           init = list(c1 = list(P0 = 17000, MSYR = .021, alpha = 3/4),
                                       c2 = list(P0 = 16000, MSYR = .017, alpha = 9/10),
                                       c3 = list(P0 = 16000, MSYR = .020, alpha = 1/4),
                                       c4 = list(P0 = 18000, MSYR = .015, alpha = 1/10)),
                           control = list(adapt_delta = .95),
                             iter = Nit, thin = 1, chains = 4)
dirichlet_sampling
stan_trace(dirichlet_sampling)
pairs(dirichlet_sampling, pars = c("P0", "MSYR", "ROI", "P1993", "alpha"))

plot(density(extract(fixed_alpha_sampling, 'P1993')$P1993),
     main = "Posterior", xlab = expression(P[1993]),
     lwd = 3
      )
lines(density(extract(dirichlet_sampling, 'P1993')$P1993), col = 2, lwd = 3)
curve(dnorm(x, m = 7800, sd = 1300), 5000, 12000, lwd = 2, lty = 2, col = 3, add = TRUE)
legend(x = "topright", legend = c("Fixed alpha", "Varying alpha", "Prior"),
       lwd = 3, lty = c(1, 1, 2), col = 1:3, bty = 'n')