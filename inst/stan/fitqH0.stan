data {
  int<lower=1> N;
  real<lower = 0, upper = 1> p[N];
  real<lower = 0, upper = 1> q[N];
  //real<lower = 0> theta[N];
  //int<lower = 0> n[N];
}
parameters {
  real<lower = 0, upper = 1> a;
  //real<lower = 0, upper = 1> b;
  real<lower = 0> sigma;
}
model {
  real phat[N];
  for (i in 1:N) phat[i] = beta_cdf(q[i], a, a);
  for (i in 1:N) {
    p[i] ~ normal(phat[i], sigma);
  }
}
