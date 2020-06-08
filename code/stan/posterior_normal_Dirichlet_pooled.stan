data {
  int<lower=0> N; // number of observations
  real Y[N]; // data
  int<lower=1> K; // number of priors to combine
  vector[K] mu; // parameters of the K priors
  vector[K] sigma_sq;
  vector<lower=0>[K] X; // hyperparameter for the alphas prior 
}
parameters {
  simplex[K] alpha;
}
transformed parameters{
  real<lower=0> vstar;
  real mustar;
  vector[K] w = alpha ./ sigma_sq;
  vstar = 1/sum(w);
  mustar =  sum(w .* mu) * vstar;
}
model {
  alpha ~ dirichlet(X);
  Y ~ normal(mustar, sqrt(vstar));
}
