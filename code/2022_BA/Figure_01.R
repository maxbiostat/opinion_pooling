source("Savchuk_data.R")
library(ggplot2)
library("logPoolR")

prior.av <- c(a0, a1, a2, a3) 
prior.bv <- c(b0, b1, b2, b3)
posterior.av <- prior.av + y
posterior.bv <- prior.bv + (n - y)
K <- length(prior.av)
expert.names <- paste0("Expert_", 0:(K-1))


theta.grid <- seq(0, 1, length.out = 1000)
prior.densities <- posterior.densities <- vector(K, mode = "list")
for(k in 1:K){
  prior.densities[[k]] <- data.frame(theta = theta.grid,
                                      dens = dbeta(theta.grid,
                                                   shape1 = prior.av[k],
                                                   shape2 = prior.bv[k]),
                                      expert = paste("expert_", k-1,
                                                     sep = ""),
                                     type = "prior")
  posterior.densities[[k]] <- data.frame(theta = theta.grid,
                                     dens = dbeta(theta.grid,
                                                  shape1 = posterior.av[k],
                                                  shape2 = posterior.bv[k]),
                                     expert = paste("expert_", k-1,
                                                    sep = ""),
                                     type = "posterior")
  
}
expert.densities.df <- rbind(
  do.call(rbind, prior.densities),
  do.call(rbind, posterior.densities)
)

expert.densities.df$type <- factor(expert.densities.df$type,
                                   levels = c("prior", "posterior"))

v.line <- data.frame(
  xintercept = y/n,
  type = "posterior"
)
v.line$type <- factor(v.line$type,
                      levels = c("prior", "posterior"))

expert_densities <- ggplot(expert.densities.df,
                           aes(x = theta, y = dens, 
                               linetype = expert, colour = expert)) + 
  geom_line(size = 2) +
  scale_linetype_manual(values = c("twodash", "dotted", "longdash", "solid"))+
  scale_colour_brewer(palette = "Paired") +
  scale_x_continuous(expression(theta), expand = c(0, 0)) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  geom_vline(data = v.line, aes(xintercept = xintercept),
  linetype = 'dashed', size = 1) +
  facet_grid(.~type) +
  theme_bw(base_size = 18) +
  theme(panel.spacing = unit(2.5, "lines"),
        legend.position = "bottom")

expert_densities
ggsave(plot = expert_densities,
       filename = "../../plots/expert_densities_Savchuk.pdf",
       width = 297, height = 210, units = "mm")

######
par_2_df <- function(parlist){
  J <- length(parlist)
  densities.list <- vector(J, mode = "list")
  for (j in 1:J){
    densities.list[[j]] <- data.frame(
      theta = theta.grid,
      dens = dbeta(theta.grid,
                   shape1 = parlist[[j]][1],
                   shape2 = parlist[[j]][2]),
      method = names(parlist)[j]
    )
  }
  densities.df <- do.call(rbind, densities.list)
  return(densities.df)
}

load("saved_data/hyperparameters_Savchuk.RData")

method.prior.densities <- par_2_df(parlist = prior.pars)
method.prior.densities$type <- "prior"
  
method.posterior.densities <- par_2_df(parlist = posterior.pars)
method.posterior.densities$type <- "posterior"
  
method.densities.df <- rbind(
  method.prior.densities,
  method.posterior.densities
)

method.densities.df$type <- factor(method.densities.df$type,
                                   levels = c("prior", "posterior"))

method_densities <- ggplot(method.densities.df,
                           aes(x = theta, y = dens, 
                               linetype = method, colour = method)) + 
  geom_line(size = 2) +
  scale_linetype_manual(values = c("twodash", "dotted",
                                   "longdash", "solid",
                                   "dotdash", "dashed")) +
  scale_colour_brewer(palette = "Set2") +
  scale_x_continuous(expression(theta), expand = c(0, 0)) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  geom_vline(data = v.line, aes(xintercept = xintercept),
             linetype = 'dashed', size = 1) +
  facet_grid(.~type) +
  theme_bw(base_size = 18) +
  theme(panel.spacing = unit(2.5, "lines"),
        legend.position = "bottom")

method_densities
ggsave(plot = method_densities,
       filename = "../../plots/method_densities_Savchuk.pdf",
       width = 297, height = 210, units = "mm")
