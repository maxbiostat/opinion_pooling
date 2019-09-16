functions{
  real L1_lpdf(real y, real df, real a, real b){
    if(is_nan(y))
      return log(0);
    else
      return student_t_lpdf((log(y + 1) - a)/b | df, 0, 1) + log(fabs(1/(b * (y + 1))));
  }
  real L2_lpdf(real y, real mu, real sd){
    if(is_nan(y))
      return log(0);
    else
      return normal_lpdf(y | mu, sd);
  }
  real qtilde_phi_lpdf(real y, real mu1, real sd1, real mu2, real sd2, real alpha){
    real ws[2] = { alpha/square(sd1), (1-alpha)/square(sd2) };
    real vstar = 1/sum(ws);
    real mstar = (ws[1]*mu1 + ws[2]*mu2)/sum(ws);
    if(is_nan(y))
      return log(0);
    else
      return normal_lpdf(y | mstar, sqrt(vstar));
      // alpha*normal_lpdf(y | mu1, sd1) + (1-alpha)*normal_lpdf(y| mu2, sd2);
  }
}
data{
  /* Prior parameters */
    real<lower=0> p0_shape; // P0 prior parameters
    real<lower=0> p0_rate;
    real<lower=0> msyr_shape; // MSYR prior parameters
    real<lower=0> msyr_rate;
    real a; // student t likelihood parameters
    real b;
    real <lower=0> roi_df;
    real phi_mu; // (pre model) prior  on P1993,  q_2(\phi)
    real <lower=0> phi_sd; 
    real ind_mu; // induced prior on P1993, q_1^*()
    real<lower=0> ind_sd; 
    /* Likelihood parameters */
    real phi_lik_mu;
    real<lower=0> phi_lik_sd;
    /*Data */
    int N; // number of time points (length of catch series)
    real<lower=0> C[N]; // catch series 
    /* Dirichlet Parameters */
    real<lower=0> a_alpha;
    real<lower=0> b_alpha;
}
parameters {
  real<lower=6400> P0; // carrying capacity
  real<lower=0> MSYR;// mean s year rate
  real<lower=0, upper=1> alpha; // logarithmic pooling weight
}
transformed parameters{
  real<lower=0> P1993;
  real ROI;
  real P[N + 1];
  P[1] = P0;
  for (t in 1:N){
    P[t + 1] = P[t] - C[t] + 1.5 * MSYR * P[t] * (1-pow(P[t]/P0, 2));
  } 
  P1993 = P[N + 1];
  ROI = pow(P1993/P[N-14], 0.06666667) - 1;// P[N-14] is P1978
}
model{
/* L_2(phi) + likROI(roi) + JacLikROI(roi) + q_1(P0, MSYR) + P1993PooledPrior  */
  target += L2_lpdf(P1993 | phi_lik_mu, phi_lik_sd);
  target += L1_lpdf(ROI | roi_df, a, b);
  target += gamma_lpdf(P0-6400 | p0_shape, p0_rate);
  target += gamma_lpdf(MSYR | msyr_shape, msyr_rate);
  target += qtilde_phi_lpdf(P1993 | ind_mu, ind_sd, phi_mu, phi_sd, alpha);
  target += beta_lpdf(alpha | a_alpha, b_alpha);
}
