library(loo)
library(rstanarm)
library(ggplot2)
library(fitdistrplus)
source("pooling_aux.r")
##################

data(Kline)
d <- Kline
d$log_pop <- log(d$population)
d$contact_high <- ifelse(d$contact=="high", 1, 0)
str(d)
N <- nrow(d) ## nobs

fit10 <-
  stan_glm(
    total_tools ~ log_pop + contact_high + log_pop * contact_high,
    family = poisson(link = "log"),
    data = d,
    prior = normal(0, 1, autoscale = FALSE),
    prior_intercept = normal(0, 100, autoscale = FALSE),
    iter = 5000,
    # seed = 2030
    seed = 666
)

loo10 <- loo(fit10, save_psis = TRUE)
loo10 <- loo(fit10, k_threshold=0.7, save_psis = TRUE)
waic10 <- waic(fit10)

fit11 <- update(fit10, formula = total_tools ~ log_pop + contact_high)
fit12 <- update(fit10, formula = total_tools ~ log_pop)
loo11 <- loo(fit11, k_threshold = 0.7, save_psis = TRUE)
loo12 <- loo(fit12, k_threshold = 0.7, save_psis = TRUE)

lpd_point <- cbind(
  loo10$pointwise[, "elpd_loo"], 
  loo11$pointwise[, "elpd_loo"], 
  loo12$pointwise[, "elpd_loo"]
)

waic11 <- waic(fit11)
waic12 <- waic(fit12)

waics <- c(
  waic10$estimates["elpd_waic", 1], 
  waic11$estimates["elpd_waic", 1], 
  waic12$estimates["elpd_waic", 1]
)

waic_wts <- exp(waics) / sum(exp(waics))
pbma_wts <- pseudobma_weights(lpd_point, BB = FALSE)
pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)

fits <- list(
  fit1 = fit10,
  fit2 = fit11,
  fit3 = fit12
)
K <- length(fits)

########
#### Optimisation
## Functions
loss <- function(y_obs, mus_i, vs_i, alpha){
  ## Expectation of squared error under the pooled predictive
 pars <- pool_par_gauss(alpha, mus_i, vs_i)
 loss <- y_obs^2 - 2*y_obs*pars[1] + pars[1]^2 + pars[2]^2
   return(loss)
}
#
compute_overall_loss <- function(dt, alpha, pars){
  N <- nrow(dt)
  if(length(pars) != N) stop("list of Gaussian parameters needs to be of same size as number of obs")
  Ls <- rep(NA, N)
  for(i in 1:N){
    Ls[i] <- loss(y_obs = dt$y[i], mus_i = pars[[i]]$mu, vs_i = pars[[i]]$v, alpha = alpha)
  }
  # cat(Ls, "\n")
  return(sum(Ls))
}
#
loss_alpha_unconstrained <- function(alpha_unc){
  compute_overall_loss(dt = d, alpha = alpha_01(alpha_unc), pars = pars.list)
}
#
parse_pars <- function(dist.fit){
  K <- length(dist.fit)
  mus <- vs <- rep(NA, K) 
  for(k in 1:K){
    mus[k] <- dist.fit[[k]]$estimate[1]
    vs[k] <- dist.fit[[k]]$estimate[2]^2
  }
  return(list(
    mu = mus,
    v = vs
  ))
}
##################
#### Getting parameters
full.postpred <- lapply(fits,
                        posterior_predict)

pars.list <- vector(N, mode = "list")
for(i in 1:N){
  pars.list[[i]] <- parse_pars(
    lapply(full.postpred, function(pp) fitdist(pp[, i], distr = dnorm, method = "mle"))
  )
}

### Optimising 
d$y <- d$total_tools

M <- 10000
quad.many.startingPoints <- matrix(rnorm(n = (K-1)*M, mean = 0, sd = 1E2), ncol = K-1, nrow = M)
many.quads <- lapply(1:M, function(i) {
  optim(quad.many.startingPoints[i, ], loss_alpha_unconstrained)
})
optimised.quads <- unlist(lapply(many.quads, function(x) x$value))

hist(optimised.quads)
abline(v = optimised.quads[which.min(optimised.quads)], lty = 2, lwd = 2)

min_quad_loss_wts <- alpha_01(many.quads[[which.min(optimised.quads)]]$par)

round(cbind(waic_wts, pbma_wts, pbma_BB_wts, stacking_wts, min_quad_loss_wts), 2)

### Pooling and plotting for each observation

for(pos in 1:nrow(d)){
  
  new.data <- d[pos, ] 
  
  mus <- pars.list[[pos]]$mu
  vs <- pars.list[[pos]]$v
  
  pars.EqualWeights <- pool_par_gauss(rep(1/K, K), mus, vs)
  pars.WAIC <- pool_par_gauss(waic_wts, mus, vs)
  pars.pBMA <- pool_par_gauss(pbma_BB_wts, mus, vs)
  pars.stacking <- pool_par_gauss(stacking_wts, mus, vs)
  pars.minQuadPool <- pool_par_gauss(min_quad_loss_wts, mus, vs)
  
  pred.dfs <- vector(length(fits), mode = "list")
  for(k in 1:K){
    pred.dfs[[k]] <- data.frame(
      y.pred.new = as.vector(posterior_predict(fits[[k]], newdata = new.data)),
      model = paste("model_", k, sep = "")
    )
  }
  all.preds <- do.call(rbind, pred.dfs)
  
  pplot <- ggplot(all.preds, aes(x = y.pred.new, colour = model, fill = model)) + 
    geom_density(alpha = .4) +
    scale_x_continuous("", expand = c(0, 0)) +
    scale_y_continuous("Density", expand = c(0, 0)) +
    geom_vline(xintercept = new.data$total_tools, linetype = "dotted") + 
    stat_function(fun = dnorm, args = list(mean = pars.minQuadPool[1],
                                           sd = pars.minQuadPool[2]), inherit.aes = FALSE, linetype = "longdash") + 
    stat_function(fun = dnorm, args = list(mean = pars.stacking[1],
                                           sd = pars.stacking[2]), inherit.aes = FALSE) + 
    theme_bw(base_size = 16) +
    ggtitle(paste("Data point", pos))
  print(pplot)
}
