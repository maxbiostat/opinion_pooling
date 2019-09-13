bowhead.SpIR <- rbind(
  read.csv("../data/output/bowhead_results_SpIR_fixed.csv"),
  read.csv("../data/output/bowhead_results_SpIR_varying.csv")
)
bowhead.MCMC <- read.csv("../data/output/bowhead_results_mcmc.csv")

bowhead.posteriors <- rbind(bowhead.SpIR, bowhead.MCMC)

## means and CIs
aggregate(bowhead.posteriors$alpha,
          by = list(bowhead.posteriors$alpha_status, bowhead.posteriors$method),
          FUN = function(x) c(mean(x, na.rm = TRUE), quantile(x, probs = c(.025, .5, .975), na.rm = TRUE)))

## Pr(alpha < 1/2)
aggregate(bowhead.posteriors$alpha,
          by = list(bowhead.posteriors$alpha_status, bowhead.posteriors$method),
          FUN = function(x) mean(x < 1/2))

bowhead.posteriors[bowhead.posteriors$alpha_status == "fixed",]$alpha <- NA

library(reshape2)

forplot.bowhead <- melt(data = bowhead.posteriors,
                        by = list(alpha_status = bowhead.posteriors$alpha_status, method = bowhead.posteriors$method),
                        variable.name = "parameter")
head(forplot.bowhead)

library(ggplot2)


posteriors_bowhead <- ggplot(forplot.bowhead, aes(x = value, colour = alpha_status, fill = alpha_status)) +
  geom_density(alpha = .4) +
  scale_x_continuous("", expand = c(0, 0)) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  facet_wrap(parameter~method, scales = "free", drop = TRUE) +
  theme_bw(base_size = 16)

posteriors_bowhead  
ggsave(posteriors_bowhead, file = "../plots/bowhead_posteriors.pdf")
