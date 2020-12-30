get_parameter_vectors_mean_cv_normal <- function(cv_correct, cv_wrong = .1){
  mv <- c(1, 2, 3, 4, 5)
  cvs <- rep(.1, length(mv))
  cvs[3] <- cv_correct
  v0s <-(cvs*mv)^2
  return(list(
    mv = as.numeric(mv),
    vv = as.numeric(v0s)
  ))
}
source("../beta_elicitator.r")
get_parameter_vectors_mean_cv_beta <- function(cv_correct, cv_wrong = .1){
  pars_0 <- elicit_beta_mean_cv(m0 = .1, cv = cv_wrong)
  pars_1 <- elicit_beta_mean_cv(m0 = .2, cv = cv_wrong)
  pars_2 <- elicit_beta_mean_cv(m0 = .5, cv = cv_correct) # big shot
  pars_3 <- elicit_beta_mean_cv(m0 = .8, cv = cv_wrong)
  pars_4 <- elicit_beta_mean_cv(m0 = .9, cv = cv_wrong)
  
  av <- unlist( c(pars_0[1], pars_1[1], pars_2[1], pars_3[1], pars_4[1]) )
  bv <- unlist( c(pars_0[2], pars_1[2], pars_2[2], pars_3[2], pars_4[2]) )
  return(list(
    av = as.numeric(av),
    bv = as.numeric(bv)
  ))
}
##############################
source("../pooling_aux.r")
##############################
as <- .001
alpha.small <- rep( (1-as)/4, 5)
alpha.small[3] <- as
# sum(alpha.small)
alpha.uniform <- rep(1/5, 5) 
#  
beta.pars.equal <- get_parameter_vectors_mean_cv_beta(cv_correct = .1, cv_wrong = .1)
beta.pars.low  <- get_parameter_vectors_mean_cv_beta(cv_correct = .0001, cv_wrong = .1)
#
normal.pars.equal <- get_parameter_vectors_mean_cv_normal(cv_correct = .1, cv_wrong = .1)
normal.pars.low  <- get_parameter_vectors_mean_cv_normal(cv_correct = .0001, cv_wrong = .1)

normal.pars.equal
normal.pars.low

pool_par(alpha = alpha.small, a = beta.pars.equal$av, b = beta.pars.equal$bv)
pool_par(alpha = alpha.small, a = beta.pars.low$av, b = beta.pars.low$bv)

pool_par(alpha = alpha.uniform, a = beta.pars.equal$av, b = beta.pars.equal$bv)
pool_par(alpha = alpha.uniform, a = beta.pars.low$av, b = beta.pars.low$bv)

#######
pool_par_gauss(alpha = alpha.small, m = normal.pars.equal$mv, v = normal.pars.equal$vv)
pool_par_gauss(alpha = alpha.small, m = normal.pars.low$mv, v = normal.pars.low$vv)

pool_par_gauss(alpha = alpha.uniform, m = normal.pars.equal$mv, v = normal.pars.equal$vv)
pool_par_gauss(alpha = alpha.uniform, m = normal.pars.low$mv, v = normal.pars.low$vv)

############
beta_expectation_alpha <- function(aright, cv = .001){
  alpha <- rep( (1-aright)/4, 5)
  alpha[3] <- aright
  pars <- get_parameter_vectors_mean_cv_beta(cv_correct = cv, cv_wrong = .1)
  par.star <- pool_par(alpha = alpha, a = pars$av, b = pars$bv)
  return(beta_mean(par.star[1], par.star[2]))
}
beta_expectation_alpha <- Vectorize(beta_expectation_alpha)
curve(beta_expectation_alpha(x), 0, .1, lwd = 3, ylab = expression(E[theta]), xlab = expression(alpha[3]) )
abline(h = 1/2, lwd = 2, lty = 2)
#########
normal_expectation_alpha <- function(aright, cv = .001){
  alpha <- rep( (1-aright)/4, 5)
  alpha[3] <- aright
  pars <- get_parameter_vectors_mean_cv_normal(cv_correct = cv, cv_wrong = .1)
  par.star <- pool_par_gauss(alpha = alpha, m = pars$mv, v = pars$vv)
  return(par.star[1])
}
normal_expectation_alpha <- Vectorize(normal_expectation_alpha)
curve(normal_expectation_alpha(x), 0, .1, lwd = 3, ylab = expression(E[theta]), xlab = expression(alpha[3]) )
abline(h = 3, lwd = 2, lty = 2)
