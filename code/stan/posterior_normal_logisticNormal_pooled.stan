data {
  int<lower=0> N; // number of trials
  int Y; // successes
  int<lower=1> K; // number of priors to combine
  vector[K] means; // mean hyperparameters for the logistic Normal (alphas prior)
  cov_matrix[K] Sigma; // variance-covariance matrix for the logistic Normal (alphas prior)
  real<lower=0> a[K]; // parameters of the K priors
  real<lower=0> b[K];
}
parameters {
real<lower=0,upper=1> theta;
vector[K] m;
}
model {
  vector[K] alpha;
  real astaraux[K];
  real bstaraux[K];
  real astar;
  real bstar; 
    for (j in 1:(K-1)) alpha[j] = exp(m[j])/(1 + sum(exp(m[1:(K-1)]))); 
    alpha[K] = 1/(1 + sum(exp(m[1:(K-1)]))); 
    for (k in 1:K){
       astaraux[k] = alpha[k]*a[k];
       bstaraux[k] = alpha[k]*b[k];
    }
  astar = sum(astaraux);
  bstar = sum(bstaraux);
  m ~ multi_normal(means, Sigma);
  theta ~ beta(astar, bstar);
  Y ~ binomial(N, theta);
}
