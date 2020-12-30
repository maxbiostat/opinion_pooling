source("pooling_aux.r")
source("beta_elicitator.r")
source("one_right_many_wrong_aux.r")

###########################

x <- 50
n <- 100
cv <- .2

expert.pars <- get_parameter_vectors_mean_cv(cv)

lapply(expert.pars, round , 2)

K <- length(expert.pars$av)

marginal.likelihoods <- rep(NA, K)
for (k in 1:K){ marginal.likelihoods[k] <- ml_beta(yi = x, ni = n, a = expert.pars$av[k], b = expert.pars$bv[k]) }

marginal.likelihoods

get_ratio(marginal.likelihoods)

alphapp <- marginal.likelihoods/sum(marginal.likelihoods)

round(alphapp, 3)

pars.pool.pp <- pool_par(alpha = alphapp, a = expert.pars$av, b = expert.pars$bv)
pars.pool.pp

beta_mean(pars.pool.pp[1], pars.pool.pp[2])
beta_mode(pars.pool.pp[1], pars.pool.pp[2])
beta_sd(pars.pool.pp[1], pars.pool.pp[2])

plot_densities(expert.pars)

curve(dbeta(x, shape1 = pars.pool.pp[1], shape2 = pars.pool.pp[2] ), lwd = 3, lty = 2, add = TRUE, col = "grey50")

