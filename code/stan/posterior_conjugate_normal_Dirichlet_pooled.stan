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
  real theta;
  real<lower=0> tau; 
}
transformed parameters{
  real<lower=0> vstar;
  real mustar;
  vector[K] w = alpha ./ sigma_sq;
  vstar = 1/sum(w);
  mustar =  sum(w .* mu) * vstar;
}
model {
  tau ~ normal(0, 5);
  alpha ~ dirichlet(X);
  theta ~ normal(mustar, sqrt(vstar));
  Y ~ normal(theta, tau);
}
