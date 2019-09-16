data {
  int<lower=0> N; // number of trials
  int Y; // successes
  int<lower=1> K; // number of priors to combine
  vector[K] means; // mean hyperparameters for the logistic Normal (alphas prior)
  cov_matrix[K] Sigma; // variance-covariance matrix for the logistic Normal (alphas prior)
  vector<lower=0>[K] a; // parameters of the K priors
  vector<lower=0>[K] b;
}
parameters {
  real<lower=0,upper=1> theta;
  vector[K] m;
}
transformed parameters{
  row_vector[K] alpha;
  real astar;
  real bstar;
  for (j in 1:(K-1)) alpha[j] = exp(m[j])/(1 + sum(exp(m[1:(K-1)]))); 
  alpha[K] = 1/(1 + sum(exp(m[1:(K-1)]))); 
  astar = alpha * a;
  bstar = alpha * b;
}
model {
  m ~ multi_normal(means, Sigma);
  theta ~ beta(astar, bstar);
  Y ~ binomial(N, theta);
}
