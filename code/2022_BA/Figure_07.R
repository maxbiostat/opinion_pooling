library(ggplot2)
library(ggthemes)
library(rstan)
library(outbreaks)

data(influenza_england_1978_school)
Ndata <- 763
sol <- influenza_england_1978_school
sol$time <- as.numeric(sol$date-min(sol$date)) + 2
sol$I <- sol$in_bed
forfit.sol <- sol
noisy_I <- forfit.sol$I/Ndata
iniTime <- 0
iniI <- 1/Ndata

devtools::source_url("https://raw.githubusercontent.com/maxbiostat/CODE/master/R/DISTRIBUTIONS/lognormal_parameter_transformation.r")

muR <- LogMean(realMean = 1.5, realSD = .25)
sdR <- sqrt(LogVar(realMean = 1.5, realSD = .25))

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

load( file = "saved_data/boarding_school_results.RData")

prediction.bands.SIR <- do.call(rbind,
                                list(predicted.varying,
                                     predicted.fixedHalf,
                                     predicted.fixedOne))

predictions_SIR <- ggplot(data = prediction.bands.SIR,
                          aes(x = time, y = post_mean,
                              colour = alpha, fill = alpha)) +
  geom_line() +
  geom_point(data = data.frame(time = epi.data$ts, I = epi.data$y),
             aes(x = time, y = I), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  scale_x_continuous("Time", expand = c(0, 0)) + 
  scale_y_continuous(expression(I(t)), expand = c(0, 0)) +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  theme_bw(base_size = 16)

predictions_SIR
ggsave(plot = predictions_SIR,
       file = "../../plots/SIR_bands.pdf")

alpha.posterior <- data.frame(alpha = extract(SIR.posterior.varying,
                                              'alpha')$alpha)

alpha_post_plot <- ggplot(data = alpha.posterior, aes (x = alpha)) +
  geom_density(alpha = .4, fill = "grey") + 
  stat_function(fun = function(x) dbeta(x, epi.data$a_alpha, epi.data$b_alpha),
                linetype = "twodash", size = 1.10) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(alpha), expand = c(0, 0)) +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  theme_bw(base_size = 20)
alpha_post_plot
ggsave(plot = alpha_post_plot,
       filename = "../../plots/posterior_alpha_SIR_example.pdf")

R0.varying <- data.frame(R0 = extract(SIR.posterior.varying, 'R0')$R0,
                         alpha = "varying") 
R0.fixedHalf <- data.frame(R0 = extract(SIR.posterior.fixedHalf, 'R0')$R0,
                           alpha = "0.5")
R0.fixedOne <- data.frame(R0 = extract(SIR.posterior.fixedOne, 'R0')$R0,
                          alpha = "1")

R0.posteriors <- do.call(rbind, list(R0.varying, R0.fixedHalf, R0.fixedOne))

R0_posterior <- ggplot(data = R0.posteriors, aes(x = R0, colour = alpha, fill = alpha)) +
  geom_density(alpha = .4) +
  geom_vline(xintercept = 0.00218/.5 * Ndata, linetype = "dotted", size = 1.01) + 
  stat_function(fun = function(x) dlnorm(x, meanlog = muR, sdlog = sdR),
                inherit.aes = FALSE, linetype = "twodash", size = 1.10) +
  stat_function(fun = function(x) dlnorm(x, meanlog = 0,
                                         sdlog = sqrt(epi.data$sigma_beta^2 + epi.data$sigma_gamma^2)),
                inherit.aes = FALSE, linetype = "F1", size = 1.10) +
  scale_x_continuous(expression(R[0]), expand = c(0, 0)) + 
  scale_y_continuous("Density", expand = c(0, 0)) + 
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  theme_bw(base_size = 16)
R0_posterior
ggsave(plot = R0_posterior,
       file = "../../plots/R0_posteriors_SIR_boardingSchool.pdf")
