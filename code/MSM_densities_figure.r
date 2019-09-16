
forplot.df.prior <- rbind(
  read.csv("../data/output/MSM_Beta_expert_densities.csv"),
  read.csv("../data/output/MSM_Gaussian_expert_densities.csv")
)

distribution_priors <- ggplot(forplot.df.prior, aes(x = phi, y = dens,
                                                  linetype = study, colour = study)) + 
  geom_line(size = 2) +
  scale_x_continuous(expression(phi), expand = c(0, 0), limits = c(0, .4)) +
  scale_y_continuous(expression(f[i](phi)), expand = c(0, 0)) +
  facet_grid(distribution~.) +
  scale_colour_brewer(palette = "Spectral") +
  theme_bw(base_size = 16)

distribution_priors
ggsave(distribution_priors, filename = "../plots/prior_densities_MSM.pdf")


distribution_priors2 <- ggplot(forplot.df.prior, aes(x = phi, y = dens,
                                                    linetype = distribution, colour = distribution)) + 
  geom_line(size = 2) +
  scale_x_continuous(expression(phi), expand = c(0, 0), limits = c(0, .4)) +
  scale_y_continuous(expression(f[i](phi)), expand = c(0, 0)) +
  facet_wrap(study~., scales = "free") +
  theme_bw(base_size = 16)

distribution_priors2
ggsave(distribution_priors2, filename = "../plots/prior_densities_byStudy_MSM.pdf")

meta <- read.csv("../data/meta_analysis_Malta_2010.csv")


forplot.df.posterior <- rbind(
  read.csv("../data/output/MSM_Beta_densities.csv"),
  read.csv("../data/output/MSM_Gaussian_densities.csv")
)

distribution_posteriors <- ggplot(forplot.df.posterior, aes(x = phi, y = dens,
                                                  linetype = method, colour = method)) + 
  geom_line(size = 2) +
  scale_x_continuous(expression(phi), expand = c(0, 0), limits = c(0, .4)) +
  scale_y_continuous(expression(pi(phi*"|"*alpha) ), expand = c(0, 0)) +
  geom_vline(xintercept = sum(meta$HIV)/sum(meta$SampleSize), linetype = "dashed", size = 1.2) +
  facet_grid(distribution~.) +
  theme_bw(base_size = 16)

distribution_posteriors
ggsave(distribution_posteriors, filename = "../plots/posterior_densities_MSM.pdf")