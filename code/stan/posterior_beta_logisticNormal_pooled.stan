data {
  int<lower=0> N; // number of trials
  int Y; // successes
  int<lower=1> K; // number of priors to combine
  vector[K] means; // mean hyperparameters for the logistic Normal (alphas prior)
  cov_matrix[K] Sigma; // variance-covariance matrix for the logistic Normal (alphas prior)
  vector<lower=0>[K] a; // parameters of the K priors
  vector<lower=0>[K] b;
}
transformed data{
  matrix[K, K] LS = cholesky_decompose(Sigma);
}
parameters {
  real<lower=0,upper=1> theta;
  vector[K] eta;
}
transformed parameters{
  row_vector[K] alpha;
  real astar;
  real bstar;
  vector[K] m = means + LS * eta;
  for (j in 1:(K-1)) alpha[j] = exp(m[j])/(1 + sum(exp(m[1:(K-1)]))); 
  alpha[K] = 1/(1 + sum(exp(m[1:(K-1)]))); 
  astar = alpha * a;
  bstar = alpha * b;
}
model {
 /* This is the model. Below we implement it without dropping any of the "constants"
  m ~ multi_normal(means, Sigma);
  theta ~ beta(astar, bstar);
  Y ~ binomial(N, theta);
  */
  target += normal_lpdf(eta | 0, 1);
  target += beta_lpdf(theta | astar, bstar);
  target += binomial_lpmf(Y | N, theta);
}
