data {
  int<lower=0> N; // number of observations
  real Y[N]; // successes
  int<lower=1> K; // number of priors to combine
  vector[K] mu; // parameters of the K priors
  vector[K] sigma_sq;
  vector[K] X; // Dirichlet parameters
}
parameters {
  simplex[K] alpha;
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
  alpha ~ dirichlet(X);
  Y ~ normal(mustar, sqrt(vstar));
}
