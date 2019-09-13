# Luiz Max Carvalho (2019)

source("pooling_aux.r")
source("beta_elicitator.r")

z <- 2.2
set.seed(666^z)


# pars_0 <- elicit_beta(m0 = .1, cv = .25)
# pars_1 <- elicit_beta(m0 = .2, cv = .25)
# pars_2 <- elicit_beta(m0 = .8, cv = .2)
# pars_3 <- elicit_beta(m0 = .8, cv = .25)
# 
# av <- unlist( c(pars_0[1], pars_1[1], pars_2[1], pars_3[1]) )
# bv <- unlist( c(pars_0[2], pars_1[2], pars_2[2], pars_3[2]) )

K <- 10
av <- runif(K, 2, 100) # seq(1.2, 50, length.out = K)
bv <- rev(av)

cbind(av, bv)

beta_mean(av, bv)
beta_sd(av, bv)
beta_mode(av, bv)

entropies <- rep(NA, K)
for(k in 1:K) entropies[k] <- entropy_beta(av[k], bv[k])

library(fields)
ES <- entropy_surface_beta(av, bv)
image.plot(ES$as, ES$bs, ES$M,
           xlab = expression(a), ylab = expression(b), horizontal = TRUE,
           cex.lab = 1.5, cex.axis = 1.5, axis.args = list(font = 2),
           legend.cex = 1.5,
           legend.lab = expression(H[pi]), main = "Entropy Beta distribution", font = 2)



# z = 1: 0.00 0.51 0.47 0.01 [0.280644]
# z = 2: 0.00 0.75 0.00 0.25 [0.280644]
# z = 3: 0.00 0.54 0.42 0.04 [0.280644]


N <- 10000 ## could increase to, say, 10000 in order to make sure, but it's fine
ent.many.startingPoints <- matrix(rnorm(n = (K-1)*N, mean = 0, sd = 100), ncol = K-1, nrow = N)
many.ents <- lapply(1:N, function(i) {
  optim(ent.many.startingPoints[i, ], optentbeta_inv, ap = av, bp = bv)
})
optimised.ents <- unlist(lapply(many.ents, function(x) x$value))

hist(optimised.ents)
abline(v = optimised.ents[which.min(optimised.ents)], lty = 2, lwd = 2)

sol <- many.ents[[which.min(optimised.ents)]]
alphaMaxEnt.opt <- alpha_01(sol$par)

sol$value
round(alphaMaxEnt.opt, 2)

( ab.MaxEnt.star <- pool_par(alphaMaxEnt.opt, av, bv) )

entropies
entropy_beta(ab.MaxEnt.star[1], ab.MaxEnt.star[2])


