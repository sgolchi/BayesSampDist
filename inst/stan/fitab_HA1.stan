// functions {
//   real fita(real phi, real nu,
//           real theta, int n) {
//     return exp(phi*sqrt(n)^(1/nu)*theta);
//           }
//   real fitb(real delta, real tau,
//           real theta, int n) {
//     return exp(-delta*sqrt(n)*theta);
//   }
// }
data {
  int<lower=1> N;
  real<lower = 1> a[N];
  //real<lower = 0> b[N];
  real theta[N];
  real theta0;
  //real<lower = 0> sigma[N];
  int<lower = 0> n[N];
}
transformed data{
  real x[N];
  for (i in 1:N) x[i] = (theta[i] - theta0)*sqrt(n[i]);
}
parameters {
  real<lower = 0> alpha1;
  real<lower = 0> alpha2;
  //real<lower = 0> alpha3;
  //real<lower = 0> alpha4;
  //real beta1;
  //real beta2;
  //real beta3;
  real<lower = 0> sigma_a;
  //real<lower = 0> sigma_b;
}
model {
  for (i in 1:N) {
    a[i] ~ lognormal(-alpha1*x[i] + alpha2*x[i]^2, sigma_a) T[1,];
}
}
