data {
  int<lower=0> N; // number of trials
  real Y[N]; // data
  real<lower=0> tau; // fixed STANDARD DEVIATION
  int<lower=1> K; // number of priors to combine
  vector[K] means; // mean hyperparameters for the logistic Normal (alphas prior)
  cov_matrix[K] Sigma; // variance-covariance matrix for the logistic Normal (alphas prior)
  vector[K] mu; // parameters of the K priors
  row_vector[K] sigma_sq;
}
transformed data{
  matrix[K, K] LS = cholesky_decompose(Sigma);
}
parameters {
  vector[K] eta;
  real theta;
}
transformed parameters{
  row_vector[K] alpha;
  real<lower=0> vstar;
  real mustar;
  vector[K] w;
  vector[K] m = means + LS * eta;
  for (j in 1:(K-1)){
    alpha[j] = exp(m[j])/(1 + sum(exp(m[1:(K-1)]))); 
    w[j] = alpha[j] / sigma_sq[j];
  } 
  alpha[K] = 1/(1 + sum(exp(m[1:(K-1)])));
  w[K] = alpha[K] / sigma_sq[K];
  vstar = 1/sum(w);
  mustar =  sum(w .* mu) * vstar;
}
model {
  // target += normal_lpdf(tau | 0, 5);
  // target += normal_lpdf(eta | 0, 1);
  // target += normal_lpdf(Y | mustar, sqrt(vstar));
  eta ~ std_normal();
  theta ~ normal(mustar, sqrt(vstar));
  Y ~ normal(theta, tau);
}
