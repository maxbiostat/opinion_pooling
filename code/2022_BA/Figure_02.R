library(logPoolR)
source("one_right_many_wrong_beta_aux_new.r")

cv <- .2
pars <- get_parameter_vectors_mean_cv(cv_correct = cv, cv_wrong = cv)
K <- length(pars$av)
Nsims <- 1000

Xa <- rep(1, K)
Xb <- rep(1/10, K)

### Sample from the hierarchical priors
alphas.dirichlet.a <- LearnBayes::rdirichlet(n = Nsims, Xa)

alphas.logisticNormal.a <- rlogisticnorm(N = Nsims,
                                         m = digamma(Xa)-digamma(Xa[K]),
                                         Sigma = constructSigma(Xa))

hpars.1.a <- apply(alphas.dirichlet.a, 1, function(alpha){
  pool_par(alpha = alpha, a = pars$av, b = pars$bv)
}) 

hpars.2.a <- apply(alphas.logisticNormal.a, 1, function(alpha){
  pool_par(alpha = alpha, a = pars$av, b = pars$bv)
})

alphas.dirichlet.b <- LearnBayes::rdirichlet(n = Nsims, Xb)

alphas.logisticNormal.b <- rlogisticnorm(N = Nsims,
                                         m = digamma(Xb)-digamma(Xb[K]),
                                         Sigma = constructSigma(Xb))

hpars.1.b <- apply(alphas.dirichlet.b, 1, function(alpha){
  pool_par(alpha = alpha, a = pars$av, b = pars$bv)
}) 

hpars.2.b <- apply(alphas.logisticNormal.b, 1, function(alpha){
  pool_par(alpha = alpha, a = pars$av, b = pars$bv)
})

#### Plotting

pdf(file = "../../plots/induced_prior.pdf", paper = "a4")

par(mfrow = c(2, 2))

plot_densities(pars, lg = FALSE, main = "Dirichlet(1)")
for (i in 1:Nsims){
  curve(dbeta(x, hpars.1.a[1, i], hpars.1.a[2, i]),
        lwd = .2, col = "grey50", add = TRUE)
}
plot_densities(pars, lg = FALSE, add = TRUE)


plot_densities(pars, lg = FALSE, main = "Logistic normal(1)")
for (i in 1:Nsims){
  curve(dbeta(x, hpars.2.a[1, i], hpars.2.a[2, i]),
        lwd = .2, col = "grey50", add = TRUE)
}
plot_densities(pars, lg = TRUE, add = TRUE)

plot_densities(pars, lg = FALSE, main = "Dirichlet(1/10)")
for (i in 1:Nsims){
  curve(dbeta(x, hpars.1.b[1, i], hpars.1.b[2, i]),
        lwd = .2, col = "grey50", add = TRUE)
}
plot_densities(pars, lg = FALSE, add = TRUE)

plot_densities(pars, lg = FALSE, main = "Logistic normal(1/10)")
for (i in 1:Nsims){
  curve(dbeta(x, hpars.2.b[1, i], hpars.2.b[2, i]),
        lwd = .2, col = "grey50", add = TRUE)
}
plot_densities(pars, lg = FALSE, add = TRUE)
dev.off()
