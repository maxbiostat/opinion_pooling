### This script implements example from page 1250 in Poole & Raftery (2000), JASA
### Original code by Gabriel Mendes (Berkeley): http://discourse.mc-stan.org/t/bayesian-melding/3011
### Implements the (unormalised) exact target of the example, mainly to demonstrate correctness
##### Copyleft (or the one to blame): Luiz Max Carvalho (2018)
fixed_alpha <- '
functions{
real fZ_exact_lpdf(real z, real ax, real bx, real ay, real by){
// notice the lack of in-built check: ay/bx < x < by/ax
real k;
real L;
real U;
k = (bx-ax)*(by-ay);
L = max({ax, ay/z});
U = min({bx, by/z});
return(log(((U *fabs(U))- (L *fabs(L)))/(2*k))) ; 
}
real q_tilde_theta_lpdf(real z, real alpha, real mX, real MX, real mY, real MY){
return(alpha * uniform_lpdf(z |0, 5) + (1-alpha)*fZ_exact_lpdf(z | mX, MX, mY, MY));
}
}
data{
real<lower=0, upper=1> alpha;
real<lower=0> max_X;
real<lower=0, upper=max_X> min_X;
real<lower=0> max_Y;
real<lower=0, upper=max_Y> min_Y;
}
parameters {
real<lower=min_X, upper=max_X> X;
real<lower=min_Y, upper=max_Y> Y;
}
model{
target += q_tilde_theta_lpdf(Y/X | alpha, min_X, max_X, min_Y, max_Y) + uniform_lpdf( X | min_X, max_X) + uniform_lpdf(Y | min_Y, max_Y) - fZ_exact_lpdf(Y/X | min_X, max_X, min_Y, max_Y);
}
generated quantities{
real<lower=min_Y/max_X, upper=max_Y/min_X> Z;
Z = Y/X;
}
'
#####################
library(rstan)
options(mc.cores = parallel::detectCores())

fixed_alpha_run <- stan(model_code = fixed_alpha,
                        data = list(alpha = .5,
                                    min_X = 2, max_X = 4,
                                    min_Y = 6, max_Y = 9),
                        iter = 5000,
                        init = list(
                          chain1 = list(X = 3.2, Y = 7),
                          chain2 = list(X = 3.5, Y = 6.5),
                          chain3 = list(X = 3, Y = 7),
                          chain4 = list(X = 2.1, Y = 8)
                        )
)
######################
fixed_alpha_run

pairs(fixed_alpha_run)
stan_trace(fixed_alpha_run)

source("../code/pooling_aux.r")
devtools::source_url("https://raw.githubusercontent.com/maxbiostat/CODE/b8473512151b0d205fd843bc291e45eac9c52e50/R/DISTRIBUTIONS/ratio_two_uniforms.r")
dZ_exact <- function(x) dpoolnorm.positive(x = x, D = list(function(x) {dunif(x, 0, 5)},
                                                           function(x) {analytic_Z(x, ax = 2, bx = 4, ay = 6, by = 9)}),
                                           alphas = c(.5, .5)
)
dZ_approx <- function(x) dpoolnorm.positive(x = x, D = list(function(x) {dunif(x, 0, 5)},
                                                            function(x) {dgamma(x, 18.3, 7.05)}),
                                            alphas = c(.5, .5)
)

Z_samples <- extract(fixed_alpha_run, 'Z')$Z
hist(Z_samples,
     probability = TRUE, main = "Z", xlab = expression(Z))
curve(dZ_exact, 1.5, 4.5, lwd = 2, add = TRUE)
curve(dZ_approx, 1.5, 4.5, lwd = 2, col = 2, add = TRUE)
legend(x = "topright", legend = c("Exact", "Gamma approximation"), col = 1:2, lwd = 2, bty = "n")

mu <- integrate( function(x) x * dZ_exact(x), 0 , Inf)
sq <- integrate( function(x) x^2 * dZ_exact(x), 0 , Inf)

mean(Z_samples); mu$value
var(Z_samples); sq$value-mu$value^2
