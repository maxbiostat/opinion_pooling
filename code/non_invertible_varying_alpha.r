### This script implements example from page 1250 in Poole & Rafetery (2000), JASA
### Original code by Gabriel Mendes (Berkeley): http://discourse.mc-stan.org/t/bayesian-melding/3011
### I provide a couple very minor improvements:
#### (i) explicit dependence on the pooling parameter (alpha);
#### (ii) an extension where alpha is given a prior rather than fixed at 1/2;
#### (iii) a simple method-of-moments estimation function to facilitate approximating the induced distribution on Z;
### Please notice this last step could easily be done (better) outside of Stan and then plugged in.
##### Copyleft (or the one to blame): Luiz Max Carvalho (2018)
fixed_alpha <- '
functions{
real [] get_MM_ests (real [] X) {
  real x_bar;
  real var_x;
  real ans[2];
  x_bar = mean(X);
  var_x = pow(sd(X), 2);
  ans[1] = pow(x_bar, 2)/var_x; // shape
  ans[2] = x_bar/var_x; // rate
  return(ans);
 }
}
data{
real<lower=0, upper=1> alpha;
int<lower=0> M; // number of samples for method of moments
real<lower=0> max_X;
real<lower=0, upper=max_X> min_X;
real<lower=0> max_Y;
real<lower=0, upper=max_Y> min_Y;
}
transformed data{
  real<lower=0> X_rep[M];
  real<lower=0> Y_rep[M];
  real<lower=0> Z_rep[M];
  real<lower=0> parms[2];
  for(m in 1:M) X_rep[m] = uniform_rng(min_X, max_X);
  for(m in 1:M) Y_rep[m] = uniform_rng(min_Y, max_Y);
  for(m in 1:M) Z_rep[m] = Y_rep[m]/X_rep[m];  
  parms = get_MM_ests(Z_rep);
}
parameters {
real<lower=0> X;
real<lower=0> Y;
}
transformed parameters {
real<lower=0> Z;
Z = Y/X;
}
model{
X ~ uniform(min_X, max_X);
Y ~ uniform(min_Y, max_Y);
target += alpha * uniform_lpdf(Z |0,5) + (1-alpha)*gamma_lpdf(Z | parms[1], parms[2]);
}'
varying_alpha <- '
functions{
real [] get_MM_ests (real [] X) {
  real x_bar;
  real var_x;
  real ans[2];
  x_bar = mean(X);
  var_x = pow(sd(X), 2);
  ans[1] = pow(x_bar, 2)/var_x; // shape
  ans[2] = x_bar/var_x; // rate
  return(ans);
 }
}
data{
int<lower=0> M; // number of samples for method of moments
real<lower=0> max_X;
real<lower=0, upper=max_X> min_X;
real<lower=0> max_Y;
real<lower=0, upper=max_Y> min_Y;
real<lower=0> a_alpha;
real<lower=0> b_alpha;
}
transformed data{
  real<lower=0> X_rep[M];
  real<lower=0> Y_rep[M];
  real<lower=0> Z_rep[M];
  real<lower=0> parms[2];
  for(m in 1:M) X_rep[m] = uniform_rng(min_X, max_X);
  for(m in 1:M) Y_rep[m] = uniform_rng(min_Y, max_Y);
  for(m in 1:M) Z_rep[m] = Y_rep[m]/X_rep[m];  
  parms = get_MM_ests(Z_rep);
}
parameters {
real<lower=0, upper=1> alpha;
real<lower=0> X;
real<lower=0> Y;
}
transformed parameters {
real<lower=0> Z;
Z = Y/X;
}
model{
X ~ uniform(min_X, max_X);
Y ~ uniform(min_Y, max_Y);
alpha ~ beta(a_alpha, b_alpha);
target += alpha * uniform_lpdf(Z |0,5) + (1-alpha)*gamma_lpdf(Z | parms[1], parms[2]);
}'
#####################
library(rstan)
options(mc.cores = parallel::detectCores())

fixed_alpha_run <- stan(model_code = fixed_alpha,
                     data = list(alpha = .5,
                                 M = 1000,
                                 min_X = 2, max_X = 4,
                                 min_Y = 6, max_Y = 9),
                     iter = 5000,
                     init = list(
                       chain1 = list(X = 3, Y = 7),
                       chain2 = list(X = 3.5, Y = 6.5),
                       chain3 = list(X = 3, Y = 7),
                       chain4 = list(X = 2.1, Y = 8)
                     ) )

alpha_varying_run <- stan(model_code = varying_alpha,
                          data = list(M = 1000,
                                      min_X = 2, max_X = 4,
                                      min_Y = 6, max_Y = 9,
                                      a_alpha = 1, b_alpha = 1),
                          iter = 5000,
                          init = list(
                            chain1 = list(X = 3, Y = 7),
                            chain2 = list(X = 3.5, Y = 6.5),
                            chain3 = list(X = 3, Y = 7),
                            chain4 = list(X = 2.1, Y = 8)
                          ) )

######################
fixed_alpha_run
alpha_varying_run

pairs(fixed_alpha_run)
stan_trace(fixed_alpha_run)

pairs(alpha_varying_run)
stan_trace(alpha_varying_run)

hist(extract(alpha_varying_run, 'alpha')$alpha,
     probability = TRUE, main = "Pooling parameter", xlab = expression(alpha))
curve(dbeta(x, 1, 1), lwd = 2, add = TRUE)