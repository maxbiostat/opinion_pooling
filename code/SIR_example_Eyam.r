library(MultiBD)
data(Eyam)
Ndata <- 350
sol <- Eyam
forfit.sol <- sol[-1, ]
noisy_I <- forfit.sol$I/Ndata + 1E-7
iniTime <- 0
iniI <- 7/Ndata

#### Elicitation
devtools::source_url("https://raw.githubusercontent.com/maxbiostat/CODE/master/R/DISTRIBUTIONS/lognormal_parameter_transformation.r")

muR <- LogMean(realMean = 2, realSD = .5)
sdR <- sqrt(LogVar(realMean = 2, realSD = .5))
qlnorm(c(.025, .975), meanlog = muR, sdlog = sdR)
curve(dlnorm(x, meanlog = muR, sdlog = sdR), 0, 5, xlab = expression(R[0]), ylab = "Density")

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
  as = 254,
  bs = 350-254,
  a_alpha = 1,
  b_alpha = 1
)
plot(epi.data$ts, epi.data$y)

#### Inference

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

SIR_code_varying <- stan_model(file = "stan/sir_simple_varying_alpha.stan")
SIR_code_fixed <- stan_model(file = "stan/sir_simple_fixed_alpha.stan")

### Varying alpha
SIR.map.varying <- optimizing(SIR_code_varying, data = epi.data, hessian = TRUE, verbose = TRUE)
SIR.posterior.varying <- sampling(SIR_code_varying, data = epi.data, chains = 4, control = list(adapt_delta = .99))
check_hmc_diagnostics(SIR.posterior.varying)
print(SIR.posterior.varying, pars = c("beta", "gamma", "S0", "R0", "sigma", "alpha"))
pairs(SIR.posterior.varying, pars = c("beta", "gamma", "S0", "R0", "sigma", "alpha"))
stan_trace(SIR.posterior.varying, pars = c("beta", "gamma", "S0", "R0", "sigma", "alpha"))
simulated_trajectories.varying <- extract(SIR.posterior.varying, 'y_rep')$y_rep
predicted.varying <- data.frame(
  time = epi.data$ts,
  lower = apply(simulated_trajectories.varying, 2, function(x) as.numeric(quantile(x, probs = .025))),
  post_mean = colMeans(simulated_trajectories.varying),
  upper = apply(simulated_trajectories.varying, 2, function(x) as.numeric(quantile(x, probs = .975))),
  alpha = "varying"
)

### Fixed alpha = 1/2
epi.data$alpha <- 1/2
SIR.map.fixedHalf <- optimizing(SIR_code_fixed, data = epi.data, hessian = TRUE, verbose = TRUE)
SIR.posterior.fixedHalf <- sampling(SIR_code_fixed, data = epi.data, chains = 4, control = list(adapt_delta = .99))
check_hmc_diagnostics(SIR.posterior.fixedHalf)
print(SIR.posterior.fixedHalf, pars = c("beta", "gamma", "S0", "R0", "sigma"))
pairs(SIR.posterior.fixedHalf, pars = c("beta", "gamma", "S0", "R0", "sigma"))
stan_trace(SIR.posterior.fixedHalf, pars = c("beta", "gamma", "S0", "R0", "sigma"))

simulated_trajectories.fixedHalf <- extract(SIR.posterior.fixedHalf, 'y_rep')$y_rep
predicted.fixedHalf <- data.frame(
  time = epi.data$ts,
  lower = apply(simulated_trajectories.fixedHalf, 2, function(x) as.numeric(quantile(x, probs = .025))),
  post_mean = colMeans(simulated_trajectories.fixedHalf),
  upper = apply(simulated_trajectories.fixedHalf, 2, function(x) as.numeric(quantile(x, probs = .975))),
  alpha = "0.5"
)

### Fixed alpha = 1
epi.data$alpha <- 1
SIR.map.fixedOne <- optimizing(SIR_code_fixed, data = epi.data, hessian = TRUE, verbose = TRUE)
SIR.posterior.fixedOne <- sampling(SIR_code_fixed, data = epi.data, chains = 4, control = list(adapt_delta = .99))
check_hmc_diagnostics(SIR.posterior.fixedOne)
print(SIR.posterior.fixedOne, pars = c("beta", "gamma", "S0", "R0", "sigma"))
pairs(SIR.posterior.fixedOne, pars = c("beta", "gamma", "S0", "R0", "sigma"))
stan_trace(SIR.posterior.fixedOne, pars = c("beta", "gamma", "S0", "R0", "sigma"))

simulated_trajectories.fixedOne <- extract(SIR.posterior.fixedOne, 'y_rep')$y_rep
predicted.fixedOne <- data.frame(
  time = epi.data$ts,
  lower = apply(simulated_trajectories.fixedOne, 2, function(x) as.numeric(quantile(x, probs = .025))),
  post_mean = colMeans(simulated_trajectories.fixedOne),
  upper = apply(simulated_trajectories.fixedOne, 2, function(x) as.numeric(quantile(x, probs = .975))),
  alpha = "1.0"
)

#### Plotting and annotating

prediction.bands.SIR <- do.call(rbind, list(predicted.varying, predicted.fixedHalf, predicted.fixedOne))

library(ggplot2)

predictions_SIR <- ggplot(data = prediction.bands.SIR, aes(x = time, y = post_mean, colour = alpha, fill = alpha)) +
  geom_line() +
  geom_point(data = data.frame(time = epi.data$ts, I = epi.data$y), aes(x = time, y = I), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  scale_x_continuous("Time", expand = c(0, 0)) + 
  scale_y_continuous(expression(I(t)), expand = c(0, 0)) + 
  theme_bw(base_size = 16)

predictions_SIR

alpha.posterior <- data.frame(alpha = extract(SIR.posterior.varying, 'alpha')$alpha)

alpha_post_plot <- ggplot(data = alpha.posterior, aes (x = alpha)) +
  geom_density(alpha = .4, fill = "grey") + 
  stat_function(fun = function(x) dbeta(x, epi.data$a_alpha, epi.data$b_alpha), linetype = "twodash", size = 1.10) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(alpha), expand = c(0, 0)) +
  theme_bw(base_size = 20)
alpha_post_plot

R0.varying <- data.frame(R0 = extract(SIR.posterior.varying, 'R0')$R0, alpha = "varying") 
R0.fixedHalf <- data.frame(R0 = extract(SIR.posterior.fixedHalf, 'R0')$R0, alpha = "0.5")
R0.fixedOne <- data.frame(R0 = extract(SIR.posterior.fixedOne, 'R0')$R0, alpha = "1")

R0.posteriors <- do.call(rbind, list(R0.varying, R0.fixedHalf, R0.fixedOne))

R0_posterior <- ggplot(data = R0.posteriors, aes(x = R0, colour = alpha, fill = alpha)) +
  geom_density(alpha = .4) +
  geom_vline(xintercept = 0.00218/.5 * Ndata, linetype = "dotted", size = 1.01) + 
  stat_function(fun = function(x) dlnorm(x, meanlog = muR, sdlog = sdR), inherit.aes = FALSE, linetype = "twodash", size = 1.10) +
  stat_function(fun = function(x) dlnorm(x, meanlog = 0, sdlog = sqrt(epi.data$sigma_beta^2 + epi.data$sigma_gamma^2)),
                inherit.aes = FALSE, linetype = "F1", size = 1.10) +
  scale_x_continuous(expression(R[0]), expand = c(0, 0)) + 
  scale_y_continuous("Density", expand = c(0, 0)) + 
  theme_bw(base_size = 16)
R0_posterior
