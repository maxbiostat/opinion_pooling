# devtools::install_github("hyunjimoon/SBC")
library(SBC)
library(LearnBayes)
library(logPoolR)
source("one_right_many_wrong_beta_aux_new.r")

cvr <- 0.2
cvw <- 0.2
pars <- get_parameter_vectors_mean_cv(cv_correct = cvr,
                                      cv_wrong = cvw)
###
Dirichlet <- 
  cmdstanr::cmdstan_model("../stan/posterior_beta_Dirichlet_pooled.stan")
backend_Dirichlet <- SBC_backend_cmdstan_sample(Dirichlet,
                                                adapt_delta = 0.99, 
                                                max_treedepth = 15)


generate_single_dataset_Dirichlet <- function(X, n) {
  K <- length(X)
  
  alphas <- as.vector( LearnBayes::rdirichlet(1, X) )  
  
  hpar <- pool_par(alpha = alphas,
                   a = pars$av,
                   b = pars$bv)
  theta <- rbeta(n = 1,
                 shape1 = hpar[1],
                 shape2 = hpar[2])
  
  y <- rbinom(n = 1, size = n, prob = theta)
  
  
  list(
    parameters = list(
      alpha = alphas,
      theta = theta,
      astar = hpar[1],
      bstar = hpar[2]
    ), generated = list(
      N = n,
      Y = y,
      K = K,
      X = X,
      a = pars$av,
      b = pars$bv
    )
  )
}

# set.seed(666)
compute <- TRUE
est_X <- 1/10
n <- 100
nrep <- 1000

if(compute){
  generator_Dirichlet <- SBC_generator_function(
    generate_single_dataset_Dirichlet, X = rep(est_X, 3), n = n
  )
  datasets_Dirichlet <- generate_datasets(generator_Dirichlet, nrep)
  save(datasets_Dirichlet,
       file =  paste0("saved_data/BetaDataSets_Dirichlet_prior_X=", est_X,
                      "_n=", n, ".RData"))
  
  results_Dirichlet <- compute_results(datasets_Dirichlet,
                                       backend_Dirichlet)
  save(results_Dirichlet,
            file = paste0("saved_data/ResultsBetaDirichlet_prior_X=",
                          est_X, "_n=", n, ".RData"))
}else{
  load(paste0("saved_data/BetaDataSets_Dirichlet_prior_X=", est_X,
                      "_n=", n, ".RData"))
  results_Dirichlet <- read.csv(
    paste0("saved_data/BetaDirichlet_prior_X=", est_X, "_n=", n, ".csv")
  )
}