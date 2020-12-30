functions{
  real pool_lpdf(real theta, vector alpha, real l1, real l2){
    real ldens = 0.0;
    ldens += alpha[1] * (normal_lpdf(theta | 2.0, 0.5) - l1);
    ldens += alpha[2] * (normal_lpdf(theta | 2.5, 2) - l2);
    ldens += alpha[3] * gamma_lpdf(theta | 16.0, 8.0);
    ldens += alpha[4] * lognormal_lpdf(theta | 0.668942, 0.7033465);
  return(ldens);
  }
}
data{
  int<lower=1> K;
  simplex[K] alpha;
}
transformed data{
  real lconst1 = normal_lccdf(0 | 2.0, 0.5);
  real lconst2 = normal_lccdf(0 | 2.5, 2.0);
}
parameters{
  real<lower=0> theta;
}
model{
  target += pool_lpdf(theta | alpha, lconst1, lconst2);
}
