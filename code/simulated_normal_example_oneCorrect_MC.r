
source("pooling_aux.r")
source("one_right_many_wrong_normal_aux.r")

require("LearnBayes")

M <- 10000

logkernel <- function(alpha, pars, data){
  y <- data$y
  n <- length(y)
  sigmasq <- data$s2
  mv <- pars$mv
  vv <- pars$vv
  ppars <- pool_par_gauss(alpha, m = mv, v = vv)
  m0 <- ppars[1]
  v0 <- ppars[2]^2
  ###
  xbar <- mean(y)
  s2 <-  sum((y-xbar)^2)
  vtilde <- n*v0 + sigmasq 
  l1 <- -s2/(2*sigmasq)
  l2 <- - (n*(xbar-m0)^2)/(2*vtilde)
  l3 <- -0.5 * ( log(vtilde) + (n-1) * log(2*pi*sigmasq) + log(2*pi))
  ans <- l1 + l2 + l3
  return(ans)
}

phi <- function(alpha, pars, data, lZ){
  y <- data$y
  n <- length(y)
  sigmasq <- data$s2
  mv <- pars$mv
  vv <- pars$vv
  ppars <- pool_par_gauss(alpha, m = mv, v = vv)
  m0 <- ppars[1]
  v0 <- ppars[2]^2
  ###
  xbar <- mean(y)
  s2 <-  sum((y-xbar)^2)
  vtilde <- n*v0 + sigmasq 
  l1 <- -s2/(2*sigmasq)
  l2 <- - (n*(xbar-m0)^2)/(2*vtilde)
  l3 <- -0.5 * ( log(vtilde) + (n-1) * log(2*pi*sigmasq) + log(2*pi))
  ans <- log(alpha) + l1 + l2 + l3 - lZ
  ans <- exp(ans)
  return(ans)
}

do_experiment <- function(cv, y, sigma, n){
  
  evidence.info <-  paste("n=", n, sep = "")
  
  pars <- get_parameter_vectors_mean_cv(cv_correct = cv)
  K <- length(pars$mv)
  X <- rep(1, K) ## Dirichlet parameter
  mock <- rep(NA, K)
  
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
  
  mpp <- pars.pool.pp[1]
  
  ## end MaL
  

  ## Dirichlet
  
  alpha.MC.dirichlet <- rdirichlet(M, X)
  
  lZ.dirichlet <- matrixStats::logSumExp(apply(alpha.MC.dirichlet, 1, logkernel, pars, data = list(y = y, s2 = fixed.sigmasq)))-log(M) 
  
  post.alpha.cred.dirichlet <- rbind(mock, mock)
  
  ## Logistic normal
  
  alpha.MC.logisticNormal <- rlogisticnorm(N = M,
                                           m = digamma(X)-digamma(X[K]),
                                           Sigma = constructSigma(X))
  lZ.logisticNormal <- matrixStats::logSumExp(apply(alpha.MC.logisticNormal, 1, logkernel, pars, data = list(y = y, s2 = fixed.sigmasq)))-log(M)
  
  post.alpha.cred.logisticNormal <-  rbind(mock, mock)
  
  ### Now "flexible stuff"
  
  Xp <- X/10

  ## Flexible Dirichlet
  
  alpha.MC.flexdirichlet <- rdirichlet(M, Xp)
  lZ.flexdirichlet <- matrixStats::logSumExp(apply(alpha.MC.flexdirichlet, 1, logkernel, pars, data = list(y = y, s2 = fixed.sigmasq)))-log(M)
  
  post.alpha.cred.flexdirichlet <-  rbind(mock, mock)
  
  ## Logistic normal 
  
  alpha.MC.flexlogisticNormal <- rlogisticnorm(N = M,
                                               m = digamma(X)-digamma(Xp[K]),
                                               Sigma = constructSigma(Xp))
  lZ.flexlogisticNormal <- matrixStats::logSumExp(apply(alpha.MC.flexlogisticNormal, 1, logkernel, pars, data = list(y = y, s2 = fixed.sigmasq)))-log(M)
  
  post.alpha.cred.flexlogisticNormal <-  rbind(mock, mock)
  ## Posterior means
  post.mean.alpha.dirichlet <-  rowMeans(apply(alpha.MC.dirichlet, 1, phi,
                                               pars = pars, data = list(y = y, s2 = fixed.sigmasq), lZ  = lZ.dirichlet))
  post.mean.alpha.flexdirichlet <-  rowMeans(apply(alpha.MC.flexdirichlet, 1, phi,
                                                   pars = pars, data = list(y = y, s2 = fixed.sigmasq), lZ  = lZ.flexdirichlet))
  post.mean.alpha.logisticNormal <- rowMeans(apply(alpha.MC.logisticNormal, 1, phi,
                                                   pars = pars, data = list(y = y, s2 = fixed.sigmasq), lZ = lZ.logisticNormal))
  post.mean.alpha.flexlogisticNormal <- rowMeans(apply(alpha.MC.flexlogisticNormal, 1, phi,
                                                       pars = pars, data = list(y = y, s2 = fixed.sigmasq), lZ = lZ.flexlogisticNormal))
  
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

set.seed(666)

truemean <- 3
fixed.sigmasq <- 1^2

ns <- c(10, 100, 500)
J <- length(ns)
dts <- lapply(ns, function(n) rnorm(n = n, mean = truemean, sd = sqrt(fixed.sigmasq)))

# cv <- .2
# y = dts[[j]]
# sigma = fixed.sigma
# n = ns[j]

per_cv <- function(cv){
  out <- lapply(seq_len(J), function(j) {
    do_experiment(cv = cv, y = dts[[j]], sigma = fixed.sigmasq, n = ns[j])
  })
  return(do.call(rbind, out))
}

cvs <- seq(0.001, 1.5, length.out = 25)
system.time(
  simu.one.correct <- do.call(rbind, lapply(cvs, per_cv))
)


library(ggplot2)

mal_ratios <- ggplot(data = subset(simu.one.correct, max_index_mal == 3), aes(x = cv, y = mal_ratio)) +
  geom_line() +
  scale_x_continuous("Coefficient of variation") + 
  scale_y_log10("Marginal likelihood ratio") + 
  facet_grid(.~data) +
  theme_bw(base_size = 16)

mal_ratios
ggsave(filename = "../plots/normal_MaLs_ratios_oneCorrect_MC.pdf", plot = mal_ratios)

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

ggsave(filename = "../plots/normal_weight_ratios_oneCorrect_MC.pdf",
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

ggsave(filename = "../plots/normal_radar_plots_oneCorrect_MC.pdf", plot = radar_alphas)

prior_mean_mal <- ggplot(data = simu.one.correct, aes(x = cv, y = mal_mean)) + ## the prior mean obtained using MaL weights
  geom_line() +
  scale_x_continuous("Coefficient of variation") + 
  scale_y_log10("Pooled prior mean") + 
  facet_grid(.~data) +
  theme_bw(base_size = 16)

prior_mean_mal
ggsave(filename = "../plots/normal_prior_mean_MaLWeights_oneCorrect_MC.pdf", plot = prior_mean_mal)