functions {
  real[] SI(real t,
  real[] y,
  real[] params,
  real[] x_r,
  int[] x_i) {
    
    real dydt[3];
    
    dydt[1] = - params[1] * y[1] * y[2];
    dydt[2] = params[1] * y[1] * y[2] - params[2] * y[2];
    dydt[3] = params[2] * y[2];
    
    return dydt;
  }
  real qtilde_phi_lpdf(real y, real mb, real sb, real mg, real sg, real mr, real sr, real alpha){
    real m1 = mb-mg;
    real v1 = square(sb) + square(sg);
    real ws[2] = { alpha/v1, (1-alpha)/square(sr) };
    real vstar = 1/sum(ws);
    real mstar = (ws[1]*m1 + ws[2]*mr) * vstar;
    return(lognormal_lpdf(y | mstar, sqrt(vstar)));
  }
}
data {
  int<lower = 0> n_obs; // Number of days sampled
  real ts[n_obs]; // Time points that were sampled
  real t0; // initial time
  real y_init; // initial measured population
  real y[n_obs]; // measured population at measurement times
  real mu_beta;
  real<lower=0> sigma_beta;
  real mu_gamma;
  real<lower=0> sigma_gamma;
  real mu_r0;
  real<lower=0> sigma_r0;
  real<lower=0> as;
  real<lower=0> bs;
  real<lower=0, upper=1> alpha;
}
transformed data{
  real x_r[2] = rep_array(0.0, 2);
  int x_i[2] = rep_array(0, 2);
}
parameters {
  real<lower=0, upper=1> r_init[2]; // initial population
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0, upper=1> S0;
  real<lower=0, upper=1> sigma;   // error scale
}
transformed parameters {
  real<lower=0> R0 = beta/gamma;
  real y_hat[n_obs, 3]; // Output from the ODE solver
  real y0[3]; // Initial conditions for both S and I
  
  y0[1] = S0;
  y0[2] = 1 - S0;
  y0[3] = 0;
  
  y_hat = integrate_ode_rk45(SI, y0, t0, ts, {beta,gamma}, x_r, x_i);
  
}
model {
  /*  beta ~ lognormal(mu_beta, sigma_beta);
  gamma ~ lognormal(mu_gamma, sigma_gamma);
  S0 ~ beta(as, bs);
  sigma ~ normal(0.5, 0.5);
  y ~ lognormal(log(y_hat[, 2]), sigma);*/
  target += beta_lpdf(S0 | as, bs);
  target += normal_lpdf(sigma | 0.5, 0.5);
  target += qtilde_phi_lpdf(R0 | mu_beta, sigma_beta, mu_gamma, sigma_beta, mu_r0, sigma_r0, alpha);
  target += lognormal_lpdf(y  | log(y_hat[, 2]), sigma);
}
generated quantities{
  real<lower=0> y_rep [n_obs];
  for (i in 1:n_obs) y_rep[i] = lognormal_rng(log(y_hat[i, 2]), sigma);
}

