data {
  int<lower=0> N; // number of observations
  real Y[N]; // successes
  int<lower=1> K; // number of priors to combine
  vector[K] mu; // parameters of the K priors
  real<lower=0> sigma_sq[K];
  real<lower=0> lambda;
}
parameters {
  real<lower=0> alpha[K];
}
transformed parameters{
  real<lower=0> vstar;
  real mustar;
  vector[K]  w;    
  for(k in 1:K) w[k] =  alpha[k] / sigma_sq[k];
  vstar = 1/sum(w);
  mustar =  sum(w .* mu) * vstar;
}
model {
  alpha ~ exponential(lambda);
  Y ~ normal(mustar, sqrt(vstar));
}
