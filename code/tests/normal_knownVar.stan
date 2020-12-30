data{
  int<lower=1> N;
  real Y[N];
  real<lower=0> tau;
  real m;
  real<lower=0> v;
}
parameters{
  real mu;
}
model{
  target += normal_lpdf(mu | m, sqrt(v));
  target += normal_lpdf(Y | mu, sqrt(tau));
}
