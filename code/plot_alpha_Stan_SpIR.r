source("pooling_aux.r")

alpha.spir <- read.csv("../data/output/alphas_spir.csv")
alpha.stan <- read.csv("../data/output/alphas_mcmc.csv")

alphas.forplot <- rbind(alpha.spir, alpha.stan)

aggregate(alphas.forplot$alpha, by = list(alphas.forplot$method), FUN = function(x) c(mean(x), quantile(x, probs = c(.025, .975))))

alpha_plot <- ggplot(alphas.forplot, aes(x = alpha, fill = method, colour = method)) +
  geom_density(alpha = .4) +
  scale_x_continuous(expression(alpha), expand = c(0, 0), breaks = number_ticks(10)) + 
  scale_y_continuous("Density", expand = c(0, 0)) +
  theme_bw(base_size = 16)

alpha_plot

alpha_plot2 <- ggplot(alphas.forplot, aes(x = method, y = alpha,  fill = method, colour = method)) +
  geom_boxplot(alpha = .4) +
  scale_x_discrete("", expand = c(0, 0)) + 
  scale_y_continuous(expression(alpha), expand = c(0, 0), breaks = number_ticks(10)) +
  theme_bw(base_size = 16)

alpha_plot2

