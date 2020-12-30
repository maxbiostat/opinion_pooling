
source("pooling_aux.r")
source("beta_elicitator.r")
source("one_right_many_wrong_beta_aux.r")

require("LearnBayes")

M <- 10000

# x <- 50
# n <- 100
# cv <- .57

logkernel <- function(alpha, pars, data){
  y <- data[1]
  n <- data[2]
  av <- pars$av
  bv <- pars$bv
  ppars <- pool_par(alpha, a = av, b = bv)
  ans <- lbeta(a = ppars[1] + y, b = ppars[2] + n-y) - lbeta(a = ppars[1], b = ppars[2])
  return(ans)
}
phi <- function(alpha, pars, data, lZ){
  y <- data[1]
  n <- data[2]
  av <- pars$av
  bv <- pars$bv
  ppars <- pool_par(alpha, a = av, b = bv)
  ans <- log(alpha) + lbeta(a = ppars[1] + y, b = ppars[2] + n-y) - lbeta(a = ppars[1], b = ppars[2]) - lZ
  ans <- exp(ans)
  return(ans)
}

do_experiment <- function(cv, x, n){
  evidence.info <-  paste("x=", x, ",", "n=", n, sep = "")
  
  pars <- get_parameter_vectors_mean_cv(cv_correct = cv)
  K <- length(pars$av)
  X <- rep(1, K) ## Dirichlet parameter
  mock <- rep(NA, K)
  
  ## Marginal likelihoods
  marginal.likelihoods <- mock
  for (k in 1:K){ marginal.likelihoods[k] <- ml_beta(yi = x, ni = n, a = pars$av[k], b = pars$bv[k]) }
  
  alphapp <- marginal.likelihoods/sum(marginal.likelihoods)
  
  pars.pool.pp <- pool_par(alpha = alphapp, a = pars$av, b = pars$bv)
  
  mpp <- beta_mean(pars.pool.pp[1], pars.pool.pp[2])
  
  ## end MaL
  
  ## Dirichlet
  
  alpha.MC.dirichlet <- rdirichlet(M, X)
  lZ.dirichlet <- matrixStats::logSumExp(apply(alpha.MC.dirichlet, 1, logkernel, pars, data = c(x, n)))-log(M)
  
  post.alpha.cred.dirichlet <- rbind(mock, mock)
  
  ## Logistic normal
  alpha.MC.logisticNormal <- rlogisticnorm(N = M,
                                           m = digamma(X)-digamma(X[K]),
                                           Sigma = constructSigma(X))
  lZ.logisticNormal <- matrixStats::logSumExp(apply(alpha.MC.logisticNormal, 1, logkernel, pars, data = c(x, n)))-log(M)
  
  post.alpha.cred.logisticNormal <-  rbind(mock, mock)
  
  ### Now "flexible stuff"
  
  Xp <- X/10
  
  ## Flexible Dirichlet
  
  alpha.MC.flexdirichlet <- rdirichlet(M, Xp)
  lZ.flexdirichlet <- matrixStats::logSumExp(apply(alpha.MC.flexdirichlet, 1, logkernel, pars, data = c(x, n)))-log(M)
  
  post.alpha.cred.flexdirichlet <-  rbind(mock, mock)
  
  ## Logistic normal 
  
  alpha.MC.flexlogisticNormal <- rlogisticnorm(N = M,
                                           m = digamma(X)-digamma(Xp[K]),
                                           Sigma = constructSigma(Xp))
  lZ.flexlogisticNormal <- matrixStats::logSumExp(apply(alpha.MC.flexlogisticNormal, 1, logkernel, pars, data = c(x, n)))-log(M)
  
  post.alpha.cred.flexlogisticNormal <-  rbind(mock, mock)
  
  ## Posterior means
  post.mean.alpha.dirichlet <-  rowMeans(apply(alpha.MC.dirichlet, 1, phi,
                                               pars = pars, data = c(x, n), lZ  = lZ.dirichlet))
  post.mean.alpha.flexdirichlet <-  rowMeans(apply(alpha.MC.flexdirichlet, 1, phi,
                                                   pars = pars, data = c(x, n), lZ  = lZ.flexdirichlet))
  post.mean.alpha.logisticNormal <- rowMeans(apply(alpha.MC.logisticNormal, 1, phi,
                                                   pars = pars, data = c(x, n), lZ = lZ.logisticNormal))
  post.mean.alpha.flexlogisticNormal <- rowMeans(apply(alpha.MC.flexlogisticNormal, 1, phi,
                                                       pars = pars, data = c(x, n), lZ = lZ.flexlogisticNormal))
  
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
  
  alpha.posteriors$mal_mean <- mpp
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

evidence <- data.frame(
  x = c(5, 50, 500, 5000),
  n = c(10, 100, 1000, 10000)
)

per_cv <- function(cv){
  out <- lapply(seq_len(nrow(evidence)), function(i) {
    do_experiment(cv = cv, x = evidence[i,]$x, n = evidence[i,]$n)
  })
  return(do.call(rbind, out))
}

cvs <- seq(0.001, .6, length.out = 25)
system.time(
  simu.one.correct <- do.call(rbind, lapply(cvs, per_cv))
)


library(ggplot2)

mal_ratios <- ggplot(data = simu.one.correct, aes(x = cv, y = mal_ratio)) +
  geom_line() +
  scale_x_continuous("Coefficient of variation") + 
  scale_y_continuous("Marginal likelihood ratio") + 
  facet_grid(.~data) +
  theme_bw(base_size = 16)

mal_ratios
ggsave(filename = "../plots/beta_MaLs_ratios_oneCorrect_MC.pdf", plot = mal_ratios)

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

ggsave(filename = "../plots/beta_weight_ratios_oneCorrect_MC.pdf",
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

ggsave(filename = "../plots/beta_radar_plots_oneCorrect_MC.pdf", plot = radar_alphas)

prior_mean_mal <- ggplot(data = simu.one.correct, aes(x = cv, y = mal_mean)) + ## the prior mean obtained using MaL weights
  geom_line() +
  scale_x_continuous("Coefficient of variation") + 
  scale_y_log10("Pooled prior mean") + 
  facet_grid(.~data) +
  theme_bw(base_size = 16)

prior_mean_mal
ggsave(filename = "../plots/beta_prior_mean_MaLWeights_oneCorrect_MC.pdf", plot = prior_mean_mal)
