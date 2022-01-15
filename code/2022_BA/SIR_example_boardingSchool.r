library(outbreaks)
library(rstan)
library(cmdstanr)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())

data(influenza_england_1978_school)
Ndata <- 763
sol <- influenza_england_1978_school
sol$time <- as.numeric(sol$date-min(sol$date)) + 2
sol$I <- sol$in_bed
forfit.sol <- sol
noisy_I <- forfit.sol$I/Ndata
iniTime <- 0
iniI <- 1/Ndata

#### Elicitation
devtools::source_url("https://raw.githubusercontent.com/maxbiostat/CODE/master/R/DISTRIBUTIONS/lognormal_parameter_transformation.r")

muR <- LogMean(realMean = 1.5, realSD = .25)
sdR <- sqrt(LogVar(realMean = 1.5, realSD = .25))
qlnorm(c(.025, .975), meanlog = muR, sdlog = sdR)
curve(dlnorm(x, meanlog = muR, sdlog = sdR), 0, 5,
      xlab = expression(R[0]), ylab = "Density")

epi.data <- list(
  n_obs = length(noisy_I),
  t0 = iniTime,
  ts = forfit.sol$time,
  y_init = iniI,
  y = noisy_I,
  mu_beta = 0,
  sigma_beta = 1,
  mu_gamma = 0,
  sigma_gamma = 1,
  mu_r0 = muR, 
  sigma_r0 = sdR,
  as = 99,
  bs = 1,
  a_alpha = 1,
  b_alpha = 1
)

#### Inference

SIR_code_varying <-
  cmdstanr::cmdstan_model(stan_file = "../stan/sir_simple_varying_alpha.stan")
SIR_code_fixed <- 
  cmdstanr::cmdstan_model(stan_file = "../stan/sir_simple_fixed_alpha.stan")

### Varying alpha
SIR.map.varying <- SIR_code_varying$optimize(data = epi.data)
SIR.posterior.varying.raw <- SIR_code_varying$sample(data = epi.data,
                                                 chains = 4,
                                                 parallel_chains = 4,
                                                 adapt_delta = .99)
SIR.posterior.varying <- stanfit(SIR.posterior.varying.raw)

check_hmc_diagnostics(SIR.posterior.varying)
print(SIR.posterior.varying,
      pars = c("beta", "gamma", "S0", "R0", "sigma", "alpha"))
pairs(SIR.posterior.varying,
      pars = c("beta", "gamma", "S0", "R0", "sigma", "alpha"))
stan_trace(SIR.posterior.varying,
           pars = c("beta", "gamma", "S0", "R0", "sigma", "alpha"))
simulated_trajectories.varying <- extract(SIR.posterior.varying, 'y_rep')$y_rep
predicted.varying <- data.frame(
  time = epi.data$ts,
  lower = apply(simulated_trajectories.varying, 2,
                function(x) as.numeric(quantile(x, probs = .025))),
  post_mean = colMeans(simulated_trajectories.varying),
  upper = apply(simulated_trajectories.varying, 2,
                function(x) as.numeric(quantile(x, probs = .975))),
  alpha = "varying"
)

### Fixed alpha = 1/2
epi.data$alpha <- 1/2
SIR.map.fixedHalf <- SIR_code_fixed$optimize(data = epi.data)
SIR.posterior.fixedHalf.raw <- SIR_code_fixed$sample(data = epi.data,
                                                   chains = 4,
                                                   parallel_chains = 4,
                                                   adapt_delta = .99)
SIR.posterior.fixedHalf <- stanfit(SIR.posterior.fixedHalf.raw)
check_hmc_diagnostics(SIR.posterior.fixedHalf)
print(SIR.posterior.fixedHalf,
      pars = c("beta", "gamma", "S0", "R0", "sigma"))
pairs(SIR.posterior.fixedHalf,
      pars = c("beta", "gamma", "S0", "R0", "sigma"))
stan_trace(SIR.posterior.fixedHalf,
           pars = c("beta", "gamma", "S0", "R0", "sigma"))

simulated_trajectories.fixedHalf <- extract(SIR.posterior.fixedHalf, 'y_rep')$y_rep
predicted.fixedHalf <- data.frame(
  time = epi.data$ts,
  lower = apply(simulated_trajectories.fixedHalf, 2,
                function(x) as.numeric(quantile(x, probs = .025))),
  post_mean = colMeans(simulated_trajectories.fixedHalf),
  upper = apply(simulated_trajectories.fixedHalf, 2,
                function(x) as.numeric(quantile(x, probs = .975))),
  alpha = "0.5"
)

### Fixed alpha = 1
epi.data$alpha <- 1
SIR.map.fixedOne <- SIR_code_fixed$optimize(data = epi.data)
SIR.posterior.fixedOne.raw <-  SIR_code_fixed$sample(data = epi.data,
                                                       chains = 4,
                                                       parallel_chains = 4,
                                                       adapt_delta = .99)
SIR.posterior.fixedOne <- stanfit(SIR.posterior.fixedOne.raw)
check_hmc_diagnostics(SIR.posterior.fixedOne)
print(SIR.posterior.fixedOne,
      pars = c("beta", "gamma", "S0", "R0", "sigma"))
pairs(SIR.posterior.fixedOne,
      pars = c("beta", "gamma", "S0", "R0", "sigma"))
stan_trace(SIR.posterior.fixedOne,
           pars = c("beta", "gamma", "S0", "R0", "sigma"))

simulated_trajectories.fixedOne <- extract(SIR.posterior.fixedOne, 'y_rep')$y_rep
predicted.fixedOne <- data.frame(
  time = epi.data$ts,
  lower = apply(simulated_trajectories.fixedOne, 2,
                function(x) as.numeric(quantile(x, probs = .025))),
  post_mean = colMeans(simulated_trajectories.fixedOne),
  upper = apply(simulated_trajectories.fixedOne, 2,
                function(x) as.numeric(quantile(x, probs = .975))),
  alpha = "1.0"
)

save(
  SIR.posterior.varying,
  SIR.posterior.fixedHalf,
  SIR.posterior.fixedOne,
  predicted.varying,
  predicted.fixedHalf,
  predicted.fixedOne,
  file = "saved_data/boarding_school_results.RData"
)

print(SIR.posterior.varying,
      pars = c("beta", "gamma", "S0", "R0", "sigma"))
print(SIR.posterior.fixedHalf,
      pars = c("beta", "gamma", "S0", "R0", "sigma"))
print(SIR.posterior.fixedOne,
      pars = c("beta", "gamma", "S0", "R0", "sigma"))
