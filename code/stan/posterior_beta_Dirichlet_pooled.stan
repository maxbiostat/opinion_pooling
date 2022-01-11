data {
  int<lower=0> N; // number of trials
  int Y; // successes
  int<lower=1> K; // number of priors to combine
  vector<lower=0>[K] X; // hyperparameter for the alphas prior 
  vector<lower=0>[K] a; // parameters of the K priors
  vector<lower=0>[K] b;
}
parameters {
real<lower=0,upper=1> theta;
simplex[K] alpha;
}
transformed parameters{
  real astar = to_row_vector(alpha) * a;
  real bstar = to_row_vector(alpha) * b;   
}
model {
/*This is the model. Below we implement it without dropping any of the "constants" */
  /*  
  alpha ~ dirichlet(X);
  theta ~ beta(astar, bstar);
  Y ~ binomial(N, theta);
  */
  target += dirichlet_lpdf(alpha | X);
  target += beta_lpdf(theta | astar, bstar);
  target += binomial_lpmf(Y | N, theta);
}
