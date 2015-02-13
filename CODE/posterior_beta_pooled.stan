data {
  int<lower=0> N; // number of trials
  int Y; // successes
  int<lower=1> K; // number of priors to combine
  vector<lower=0>[K] X; // hyperparameter for the alphas prior 
  real<lower=0> a[K]; // parameters of the K priors
  real<lower=0> b[K];
}
parameters {
real<lower=0,upper=1> theta;
simplex[K] alpha;
}
model {
  real astaraux[K];
  real bstaraux[K];
  real astar;
  real bstar; 
    for (k in 1:K){
       astaraux[k] <- alpha[k]*a[k];
       bstaraux[k] <- alpha[k]*b[k];
    }
  astar <- sum(astaraux);
  bstar <- sum(bstaraux);
  alpha ~ dirichlet(X);
  theta ~ beta(astar, bstar);
  Y ~ binomial(N, theta);
}