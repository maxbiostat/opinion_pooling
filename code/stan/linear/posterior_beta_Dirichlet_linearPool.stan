functions{
   real beta_mixture_lpdf(real theta, real[] as, real[] bs,  real[] alpha) {
    int K = size(as);
    real lmix[K];
    if(size(alpha) != K) reject("alpha is not the right dimension");
    for(k in 1:K){
      lmix[k] = log(alpha[k]) + beta_lpdf(theta | as[k], bs[k]);
    }
    return(log_sum_exp(lmix));
  }
}
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
/*This is the model. Below we implement it without dropping any of the "constants" */
  /*  
  alpha ~ dirichlet(X);
  theta ~ beta(astar, bstar);
  Y ~ binomial(N, theta);
  */
  target += dirichlet_lpdf(alpha | X);
  target += beta_mixture_lpdf(theta | a, b, to_array_1d(alpha));
  target += binomial_lpmf(Y | N, theta);
}
