source("../pooling_aux.r")
source("../beta_elicitator.r")

get_parameter_vectors <- function(cv_correct, cv_wrong = .2){
  pars_0 <- elicit_beta(m0 = .1, cv = cv_wrong)
  pars_1 <- elicit_beta(m0 = .2, cv = cv_wrong)
  pars_2 <- elicit_beta(m0 = .5, cv = cv_correct) # big shot
  pars_3 <- elicit_beta(m0 = .8, cv = cv_wrong)
  pars_4 <- elicit_beta(m0 = .9, cv = cv_wrong)
  
  av <- unlist( c(pars_0[1], pars_1[1], pars_2[1], pars_3[1], pars_4[1]) )
  bv <- unlist( c(pars_0[2], pars_1[2], pars_2[2], pars_3[2], pars_4[2]) )
  return(list(
    av = as.numeric(av),
    bv = as.numeric(bv)
  ))
}

plot_densities <- function(pars){
  av <- pars$av
  bv <- pars$bv
  K <- length(av)
  curve(dbeta(x, av[1], bv[1]), lwd = 2)
  for(k in 2:K){
    curve(dbeta(x, av[k], bv[k]), lwd = 2, col = k, add = TRUE)
  }
}

pars <- get_parameter_vectors(.15)
K <- length(pars$av)

plot_densities(pars)

opt_meanBeta_inv <- function(alpha.inv, ap, bp, target = 1/2, verbose = FALSE){
  alpha <- alpha_01(alpha.inv)
  astar <- sum(alpha*ap)
  bstar <- sum(alpha*bp)
  est <- beta_mean(astar, bstar)
  if(verbose) cat("attained mean:", est, " attained sd:", beta_sd(astar, bstar), "\n")
  return((est-target)^2)
}

( opt0 <- optim(alpha_real(rep(1/K, K)), opt_meanBeta_inv , ap = pars$av, bp = pars$bv) )

( opt1 <- optim(alpha_real(rep(1/K, K)), opt_meanBeta_inv , ap = pars$av, bp = pars$bv,
               method = "L-BFGS-B", lower = rep(-200, K) , upper = rep(200, K)) )

round(alpha_01(opt0$par), 3)
round(alpha_01(opt1$par), 3)

p3 <- .95
the.weights <- rep( (1-p3)/(K-1), K)
the.weights[3] <- p3

the.weights

sum(the.weights)
opt_meanBeta_inv(alpha.inv = alpha_real(the.weights), ap = pars$av, bp = pars$bv, verbose = TRUE)
opt_meanBeta_inv(alpha.inv = opt1$par, ap = pars$av, bp = pars$bv, verbose = TRUE)
