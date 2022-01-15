library(cmdstanr)
library(rstan)
library(Ternary)
library(logPoolR)

source("one_right_many_wrong_beta_aux_new.r")

stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())

compiled.dirichlet <-
  cmdstanr::cmdstan_model("../stan/posterior_beta_Dirichlet_pooled.stan")
compiled.logisticNormal <-
  cmdstanr::cmdstan_model("../stan/posterior_beta_logisticNormal_pooled.stan")

run_scenario <- function(sampleMean, n, cvr, cvw, est_X){
  
  hyperpars <- get_parameter_vectors_mean_cv(cv_correct = cvr,
                                             cv_wrong = cvw)    
  
  prior.means <- logPoolR:::beta_mean(hyperpars$av, hyperpars$bv)
  
  K <- length(hyperpars$av)
  
  y <-  round(n * sampleMean)
  
  marginal.likelihoods <- densities.true <- rep(NA, K)
  
  for (k in 1:K){
    marginal.likelihoods[k] <- logPoolR:::ml_beta(yi = y,
                                                  ni = n,
                                                  a = hyperpars$av[k],
                                                  b = hyperpars$bv[k])
    densities.true[k] <- dbeta(x = sampleMean,
                               shape1 = hyperpars$av[k],
                               shape = hyperpars$bv[k])
  }
  
  alphas.bma <- marginal.likelihoods/sum(marginal.likelihoods)
  
  hpar.bma <- logPoolR::pool_par(alpha = alphas.bma,
                                 a = hyperpars$av,
                                 b = hyperpars$bv)
  
  ### Fit hierarchical with Dirichlet prior
  dirichlet.data <- list(Y = y,
                         N = n,
                         X = est_X,
                         K = K,
                         a = hyperpars$av,
                         b = hyperpars$bv)
  sink('/dev/null')
  fit.dirichlet <- stanfit(compiled.dirichlet$sample(
    data = dirichlet.data,
    max_treedepth = 10,
    adapt_delta = .99,
    iter_warmup = 2000,
    iter_sampling = 2000,
    refresh = 0,
    parallel_chains = 4
  ))
  sink()
  
  posterior.alphas.Dirichlet <- 
    extract(fit.dirichlet, 'alpha')$alpha
  
  logisticNormal.data <- dirichlet.data
  logisticNormal.data$means <- digamma(est_X)-digamma(est_X[K])
  logisticNormal.data$Sigma <- constructSigma(est_X)
  
  sink('/dev/null')
  fit.logisticNormal <-
    stanfit(compiled.logisticNormal$sample(
      data = logisticNormal.data,
      max_treedepth = 10,
      adapt_delta = .99,
      iter_warmup = 2000,
      iter_sampling = 2000,
      refresh = 0,
      parallel_chains = 4)
    )
  sink()
  
  posterior.alphas.logisticNormal <- 
    extract(fit.logisticNormal, 'alpha')$alpha
  
  #### putting the output together
  
  bma <- data.frame(
    mean = c(NA, alphas.bma, hpar.bma),
    lwr = NA,
    median = NA,
    upr = NA,
    method = "BMA"
  )
  dirichlet.res <- data.frame(
    summary(fit.dirichlet,
            pars = c("theta", "alpha",
                     "astar", "bstar"))$summary[, c(1, 4, 6, 8)],
    method = "Dirichlet"
  )
  colnames(dirichlet.res) <- c("mean", "lwr", "median", "upr", "method")
  
  logisticNormal.res <- data.frame(
    summary(fit.logisticNormal,
            pars = c("theta", "alpha",
                     "astar", "bstar"))$summary[, c(1, 4, 6, 8)],
    method = "Logistic_normal"
  )
  colnames(logisticNormal.res) <- c("mean", "lwr", "median", "upr", "method")
  
  inference <- do.call(rbind, 
                       list(
                         bma,
                         dirichlet.res,
                         logisticNormal.res
                       ))
  
  out <- data.frame(
    par = rownames(dirichlet.res),
    inference,
    xbar = sampleMean,
    n = n, 
    cv_right = cvr,
    cv_wrong = cvw,
    X = est_X[1])
  rownames(out) <- NULL
  
  fName <- paste0("saved_data/oracle_beta/oracleSimulationBeta_",
                  cvr, "_", cvw,
                  "_x_bar=", sampleMean,
                  "_n=", n,
                  "_X=", est_X[1], ".RData")
  save(hyperpars,
       alphas.bma,
       hpar.bma,
       posterior.alphas.Dirichlet,
       posterior.alphas.logisticNormal,
       file = fName)
  
  return(out)
}

xbs <- 1:9/10
sizes <- c(10, 100, 1E4)

grid <- expand.grid(xbar = xbs, n = sizes)
J <- nrow(grid)

simus <- do.call(rbind,
                 lapply(1:J, function(j){
                   run_scenario(sampleMean = grid[j, 1],
                                n = grid[j, 2],
                                cvr = .2,
                                cvw = .2,
                                est_X = rep(1, 3))
                 }))

write.csv(simus,
          file = "saved_data/oracle_beta/oracle_Beta_estX=1.csv")
