data {
  int<lower=1> N;
  real<lower = 0, upper = 1> a[N];
  //real<lower = 0> b[N];
  int<lower = 0> n[N];
}
parameters {
  real<lower = 0> alpha1;
  real<lower = 0> alpha2;
  real<lower = 0> sigma_a;
}
model {
  for (i in 1:N) {
    a[i] ~ lognormal(-alpha1/n[i] + alpha2/n[i]^2 , sigma_a);
}
}
