
source("pooling_aux.r")
source("one_right_many_wrong_normal_aux.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

compiled.dirichlet <- stan_model("stan/posterior_fixedVariance_normal_Dirichlet_pooled.stan")
compiled.logisticNormal <- stan_model("stan/posterior_fixedVariance_normal_logisticNormal_pooled.stan")

do_experiment <- function(cv, y, sigma, n){
  
  evidence.info <-  paste("n=", n, sep = "")
  
  pars <- get_parameter_vectors_mean_cv(cv_correct = cv)
  K <- length(pars$mv)
  X <- rep(1, K) ## Dirichlet parameter
  
  ## Marginal likelihoods
  require(matrixStats)
  log.marginal.likelihoods <- rep(NA, K)
  for (k in 1:K){ 
    log.marginal.likelihoods[k] <- normal_mean_marg_like(y = y, sigma = sigma,
                                                         m = pars$mv[k], v = pars$vv[k],
                                                         log = TRUE)
  }
  logS <- matrixStats::logSumExp(log.marginal.likelihoods)
  alphapp <- exp(log.marginal.likelihoods-logS)
  marginal.likelihoods <- exp(log.marginal.likelihoods)
  
  pars.pool.pp <- pool_par_gauss(alpha = alphapp, m = pars$mv, v = pars$vv)
  
  
  ## end MaL
  
  gaussdata.stan <- list(Y = y,
                         tau = sigma,
                         X = X,
                         N = length(y),
                         mu = pars$mv,
                         sigma_sq = pars$vv, 
                        means = digamma(X)-digamma(X[K]),
                        Sigma = constructSigma(X),
                        nu = 1,
                        s = 1,
                        N = n, K = K, a = pars$mv, b = pars$vv)
  
  ## Dirichlet
  
  dirichlet.posterior <- sampling(compiled.dirichlet, data = gaussdata.stan,
                                  refresh = 0, control = list(adapt_delta = .99, max_treedepth = 15))
  alphas.dirichlet <- extract(dirichlet.posterior, 'alpha')$alpha
  post.alpha.cred.dirichlet <- apply(alphas.dirichlet, 2, quantile, probs = c(.025, .975))
  
  ## Logistic normal
  
  logisticNormal.posterior <- sampling(compiled.logisticNormal, data = gaussdata.stan,
                                       refresh = 0, control = list(adapt_delta = .99, max_treedepth = 15))
  alphas.logisticNormal <- extract(logisticNormal.posterior, 'alpha')$alpha
  post.alpha.cred.logisticNormal <- apply(alphas.logisticNormal, 2, quantile, probs = c(.025, .975))
  
  ### Now "flexible stuff"
  
  Xp <- X/10
  gaussdata.stan$X <- Xp
  gaussdata.stan$means <-  digamma(Xp)-digamma(Xp[K])
  gaussdata.stan$Sigma <- constructSigma(Xp)
  
  ## Flexible Dirichlet
  
  flexdirichlet.posterior <- sampling(compiled.dirichlet, data = gaussdata.stan,
                                      refresh = 0,
                                      control = list(adapt_delta = .99, max_treedepth = 15))
  alphas.flexdirichlet <- extract(flexdirichlet.posterior, 'alpha')$alpha
  post.alpha.cred.flexdirichlet <- apply(alphas.flexdirichlet, 2, quantile, probs = c(.025, .975))
  
  ## Logistic normal 
  
  flexlogisticNormal.posterior <- sampling(compiled.logisticNormal, data = gaussdata.stan,
                                           refresh = 0, control = list(adapt_delta = .99, max_treedepth = 15))
  alphas.flexlogisticNormal <- extract(flexlogisticNormal.posterior, 'alpha')$alpha
  post.alpha.cred.flexlogisticNormal <- apply(alphas.flexlogisticNormal, 2, quantile, probs = c(.025, .975))
  
  ## Posterior means
  post.mean.alpha.dirichlet <- colMeans(alphas.dirichlet)
  post.mean.alpha.flexdirichlet <- colMeans(alphas.flexdirichlet)
  post.mean.alpha.logisticNormal <- colMeans(alphas.logisticNormal)
  post.mean.alpha.flexlogisticNormal <- colMeans(alphas.flexlogisticNormal)
  
  alpha.posteriors <- data.frame(
    alpha = c(post.mean.alpha.dirichlet, post.mean.alpha.flexdirichlet,
              post.mean.alpha.logisticNormal, post.mean.alpha.flexlogisticNormal),
    lwr = c(post.alpha.cred.dirichlet[1, ], post.alpha.cred.flexdirichlet[1, ],
            post.alpha.cred.logisticNormal[1, ], post.alpha.cred.flexlogisticNormal[1, ]),
    upr = c(post.alpha.cred.dirichlet[2, ], post.alpha.cred.flexdirichlet[2, ],
            post.alpha.cred.logisticNormal[2, ], post.alpha.cred.flexlogisticNormal[2, ]),
    expert = rep(paste("expert_", 0:(K-1), sep = ""), 4),
    prior = rep(c("Dirichlet", "Flexible_Dirichlet", "Logistic-normal", "Flexible_Logistic-normal"), each = K)
  )
  alpha.posteriors$data <- evidence.info
  alpha.posteriors$cv <- cv
  
  alpha.posteriors$mal_mean <- pars.pool.pp[1]
  alpha.posteriors$mal_ratio <- get_ratio(marginal.likelihoods)
  
  alpha.posteriors$alpha_ratio <- NA
  
  alpha.posteriors[which(alpha.posteriors$prior == "Dirichlet"), ]$alpha_ratio <- get_ratio(post.mean.alpha.dirichlet)
  alpha.posteriors[which(alpha.posteriors$prior == "Flexible_Dirichlet"), ]$alpha_ratio <- get_ratio(post.mean.alpha.flexdirichlet)
  alpha.posteriors[which(alpha.posteriors$prior == "Logistic-normal"), ]$alpha_ratio <- get_ratio(post.mean.alpha.logisticNormal)
  alpha.posteriors[which(alpha.posteriors$prior == "Flexible_Logistic-normal"), ]$alpha_ratio <- get_ratio(post.mean.alpha.flexlogisticNormal)
  
  alpha.posteriors$max_index_mal <- which.max(marginal.likelihoods) 
  
  alpha.posteriors$max_index_priors <- NA
  
  alpha.posteriors[which(alpha.posteriors$prior == "Dirichlet"), ]$max_index_priors <- which.max(post.mean.alpha.dirichlet)
  alpha.posteriors[which(alpha.posteriors$prior == "Flexible_Dirichlet"), ]$max_index_priors <- which.max(post.mean.alpha.flexdirichlet)
  alpha.posteriors[which(alpha.posteriors$prior == "Logistic-normal"), ]$max_index_priors <- which.max(post.mean.alpha.logisticNormal)
  alpha.posteriors[which(alpha.posteriors$prior == "Flexible_Logistic-normal"), ]$max_index_priors <- which.max(post.mean.alpha.flexlogisticNormal)
  
  return(alpha.posteriors)  
}

set.seed(666)

truemean <- 3
fixed.sigma <- 1

ns <- c(10, 100, 500)
J <- length(ns)
dts <- lapply(ns, function(n) rnorm(n = n, mean = truemean, sd = fixed.sigma))

# cv <- .2
# y = dts[[j]]
# sigma = fixed.sigma
# n = ns[j]

per_cv <- function(cv){
  out <- lapply(seq_len(J), function(j) {
    do_experiment(cv = cv, y = dts[[j]], sigma = fixed.sigma, n = ns[j])
  })
  return(do.call(rbind, out))
}

cvs <- seq(0.001, .1, length.out = 25)
system.time(
  simu.one.correct <- do.call(rbind, lapply(cvs, per_cv))
)


library(ggplot2)

mal_ratios <- ggplot(data = simu.one.correct, aes(x = cv, y = mal_ratio)) +
  geom_line() +
  scale_x_continuous("Coefficient of variation") + 
  scale_y_log10("Marginal likelihood ratio") + 
  facet_grid(.~data) +
  theme_bw(base_size = 16)

mal_ratios
ggsave(filename = "../plots/normal_MaLs_ratios_oneCorrect.pdf", plot = mal_ratios)

weight_ratios <- ggplot(data = simu.one.correct, aes(x = cv, y = alpha_ratio, colour = prior)) +
  geom_line(size = 1.5) +
  scale_x_continuous("Coefficient of variation") + 
  scale_y_continuous("Weight ratio") + 
  facet_grid(.~data) +
  theme_bw(base_size = 16)

weight_ratios

weight_ratios_expert2 <- ggplot(data = subset(simu.one.correct, max_index_priors == 3),
                                aes(x = cv, y = alpha_ratio, colour = prior)) +
  geom_line(size = 1.5) +
  scale_x_continuous("Coefficient of variation") + 
  scale_y_continuous("Weight ratio") + 
  facet_grid(.~data) +
  theme_bw(base_size = 16)

weight_ratios_expert2

ggsave(filename = "../plots/normal_weight_ratios_oneCorrect.pdf",
       plot = weight_ratios_expert2)


radar_alphas <- ggplot(data = subset(simu.one.correct, cv %in% cvs[round(seq(1, length(cvs), length.out = 3))]),
                       aes(x = expert, y = alpha, group = prior, colour = prior, fill = prior)) +
  geom_point() +
  geom_polygon(alpha = 0.4) +
  facet_grid(cv ~ data) +
  theme_bw(base_size = 16) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                     breaks = number_ticks(10)) + 
  coord_radar() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )
radar_alphas

ggsave(filename = "../plots/normal_radar_plots_oneCorrect.pdf", plot = radar_alphas)
