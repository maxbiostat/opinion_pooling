library(rstan)
library(cmdstanr)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())


CtData <- read.csv("../../data/bowhead_kills_1848-2002.csv")
CtData <- subset(CtData, year <= 1992)
bowhead.data <- list(
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

fixed_alpha <- cmdstanr::cmdstan_model("../stan/bowhead_fixed_alpha.stan")

fixed.alpha.mcmc.raw <- fixed_alpha$sample(
  data = bowhead.data,
  chains = 4,
  parallel_chains = 4,
  adapt_delta = .99,
  max_treedepth = 15,
  init = list(c1 = list(P0 = 17000, MSYR = .021),
                          c2 = list(P0 = 16000, MSYR = .017),
                          c3 = list(P0 = 16000, MSYR = .020),
                          c4 = list(P0 = 18000, MSYR = .015))
) 
# init = ,
fixed.alpha.mcmc <- stanfit(fixed.alpha.mcmc.raw)
stan_trace(fixed.alpha.mcmc,
           pars = c("P0", "MSYR", "ROI", "P1993"))
pairs(fixed.alpha.mcmc,
      pars = c("P0", "MSYR", "ROI", "P1993"))

#####################
# Now sampling with alpha
bowhead.data2 <- bowhead.data
bowhead.data2$alpha <- NULL
bowhead.data2$a_alpha <- 1
bowhead.data2$b_alpha <- 1

bowhead_dirichlet <- cmdstanr::cmdstan_model("../stan/bowhead_varying_alpha.stan")

varying.alpha.mcmc.raw <- bowhead_dirichlet$sample(
  data = bowhead.data2,
  chains = 4,
  parallel_chains = 4,
  adapt_delta = .99,
  max_treedepth = 15,
  init = list(c1 = list(P0 = 17000, MSYR = .021),
              c2 = list(P0 = 16000, MSYR = .017),
              c3 = list(P0 = 16000, MSYR = .020),
              c4 = list(P0 = 18000, MSYR = .015))
)
varying.alpha.mcmc <- stanfit(varying.alpha.mcmc.raw)
stan_trace(varying.alpha.mcmc)
pairs(varying.alpha.mcmc,
      pars = c("P0", "MSYR", "ROI", "P1993", "alpha"))

plot(density(extract(fixed.alpha.mcmc, 'P1993')$P1993),
     main = "Posterior", xlab = expression(P[1993]),
     lwd = 3)
lines(density(extract(varying.alpha.mcmc, 'P1993')$P1993), col = 2, lwd = 3)
curve(dnorm(x, m = 7800, sd = 1300), 5000, 12000, lwd = 2,
      lty = 2, col = 3, add = TRUE)
legend(x = "topright",
       legend = c("Fixed alpha", "Varying alpha", "Prior"),
       lwd = 3, lty = c(1, 1, 2), col = 1:3, bty = 'n')


### Exporting stuff
df.fixed <- data.frame(
  P0 = extract(fixed.alpha.mcmc, 'P0')$P0,
  MSYR = extract(fixed.alpha.mcmc, 'MSYR')$MSYR,
  ROI = extract(fixed.alpha.mcmc, 'ROI')$ROI,
  P1993 = extract(fixed.alpha.mcmc, 'P1993')$P1993,
  alpha = 1/2,
  alpha_status = "fixed",
  method = "MCMC"
)

df.varying <- data.frame(
  P0 = extract(varying.alpha.mcmc, 'P0')$P0,
  MSYR = extract(varying.alpha.mcmc, 'MSYR')$MSYR,
  ROI = extract(varying.alpha.mcmc, 'ROI')$ROI,
  P1993 = extract(varying.alpha.mcmc, 'P1993')$P1993,
  alpha = extract(varying.alpha.mcmc, 'alpha')$alpha,
  alpha_status = "varying",
  method = "MCMC"
)

Theta.posteriors <- rbind(
  df.fixed,
  df.varying
)

write.csv(Theta.posteriors, 
          file = "saved_data/bowhead_results_mcmc.csv",
          row.names = FALSE)
