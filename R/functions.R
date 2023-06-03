#' Quantile matching under the null
#'
#' @export
#' @param ps Numeric vector of posterior probabilities (empirical distribution of the test-statistic).
#' @param q Numeric vector of tail probabilities corresponsing to the quantiles to be matched.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return A sample of shape parameter values for a beta fit with equal shape and scale parameters
#'
qmatchH0 = function(ps, q = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.98, 0.985, 0.99), ...) {
  p = p0 = NULL
  for (i in 1:length(q)) {
    p[i] = mean(ps<q[i])
  }
  i1 = ifelse(length(which(p>0))>0, min(which(p>0)), 1)
  i2 = ifelse(length(which(p>1))>0, max(which(p>1)), length(q))
  q1 = q[i1:length(q)]
  p = p[i1:i2]
  standata = list(N = length(p), p = p, q = q1)
  fit0 = rstan::sampling(stanmodels$fitqH0, data = standata, chains = 2, refresh = 0, ...)
  fitss = rstan::extract(fit0)
  a = sample(fitss$a, 100, replace = T)
  return(a)
}

#' Quantile matching under the alternative
#'
#' @export
#' @param ps Numeric vector of posterior probabilities (empirical distribution of the test-statistic).
#' @param q Numeric vector of tail probabilities corresponsing to the quantiles to be matched.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return A sample of shape parameter values for a beta fit with shape=1/scale parameters
#'
qmatchHA = function(ps, q = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.98, 0.985, 0.99), ...) {
  p = p0 = NULL
  for (i in 1:length(q)) {
    p[i] = mean(ps<q[i])
  }
  i1 = ifelse(length(which(p>0))>0, min(which(p>0)), 1)
  i2 = ifelse(length(which(p>1))>0, max(which(p>1)), length(q))
  q1 = q[i1:length(q)]
  p = p[i1:i2]
  standata = list(N = length(p), p = p, q = q1)
  fit0 = rstan::sampling(stanmodels$fitqHA, data = standata, chains = 2, refresh = 0, ...)
  fitss = rstan::extract(fit0)
  a = sample(fitss$a, 100, replace = T)
  return(a)
}

#' Beta parameters fit under the null
#'
#' @export
#' @param dfH0 a dataframe with two columns 'a' and 'n', containing the 'a' samples from qmatchH0 with the corresponding n.
#' @return list of estimated parameters (posterior means)
#'
abfith0 = function(dfH0) {
  standata = list(N = nrow(dfH0), a = dfH0$a, n = dfH0$n)
  fit = rstan::sampling(stanmodels$fitab_H0, data = standata, chains = 1)
  fitss = rstan::extract(fit)
  alpha1_hat = mean(fitss$alpha1)
  alpha2_hat = mean(fitss$alpha2)
  siga_hat = mean(fitss$sigma_a)
  out = list(alpha1 = alpha1_hat, alpha2 = alpha2_hat, sig0 = siga_hat)
  return(out)
}

#' Beta parameters fit under the alternative
#'
#' @export
#' @param dfHA a dataframe with two columns 'a', 'n', and 'theta' containing the 'a' samples from qmatchHA with the corresponding n and theta>0.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return list of estimated parameters (posterior means)
#'
abfitha = function(dfHA) {
  standata = list(N = nrow(dfHA), a = dfHA$a, theta0 = 0, theta = dfHA$theta, n = dfHA$n)
  fit = rstan::sampling(stanmodels$fitab_HA1, data = standata, chains = 1)
  fitss = rstan::extract(fit)
  alpha1_hat = mean(fitss$alpha1)
  alpha2_hat = mean(fitss$alpha2)
  siga_hat = mean(fitss$sigma_a)
  out = list(beta1 = alpha1_hat, beta2 = alpha2_hat, sig1 = siga_hat)
  return(out)
}

#' t1e estimate
#'
#' @export
#' @param n sample size.
#' @param u decision threshold.
#' @param pars list of beta parameter estimates under the null (alpha1, alpha2, sig0)
#' @return vector of posterior mean and 95% CIs for t1e
#'
pEst0 = function(n, u, pars) {
  asamp = exp(rnorm(1000, -pars$alpha1/(n) + pars$alpha2/n^2, pars$sig0))
  ps_hat = NULL
  for (i in 1:length(asamp)) {
    ps_hat = c(ps_hat, 1 - pbeta(u, asamp[i], asamp[i]))
  }
  out = c(psh = mean(ps_hat), psl = quantile(ps_hat, 0.025), psu = quantile(ps_hat, 0.975))
  return(out)
}

#' power estimate
#'
#' @export
#' @param n sample size.
#' @param theta effect size (linear scale).
#' @param u decision threshold.
#' @param pars list of beta parameter estimates under the alternative (beta1, beta2, sig1)
#' @return vector of posterior mean and 95% CIs for power
#'
pEstA = function(n, theta, u, pars) {
  theta  = -abs(theta)
  x = (theta)*sqrt(n)
  asamp = exp(rnorm(1000, -pars$beta1*x + pars$beta2*x^2, pars$sig1))
  ps_hat = NULL
  for (i in 1:length(asamp)) {
    ps_hat = c(ps_hat, 1 - pbeta(u, asamp[i], 1/asamp[i]))
  }
  out = c(psh = mean(ps_hat), psl = quantile(ps_hat, 0.025), psu = quantile(ps_hat, 0.975))
  return(out)
}

#' t1e pred
#'
#' @export
#' @param testcases test/prediction set: a vector of sample sizes.
#' @param u decision threshold.
#' @param pars list of beta parameter estimates under the null (alpha1, alpha2, sig0)
#' @return data frame of posterior mean and 95% CIs for t1e in test set
#'
pred0 = function(testcases, u, pars) {
  psh = sapply(testcases, pEst0, u = u, pars = pars)
  df = as.data.frame(t(simplify2array(psh)))
  df = data.frame(cbind(df, n = testcases))
  names(df) = c('psh', 'psl', 'psu', 'n')
  return(df)
}

#' power pred
#'
#' @export
#' @param testcases test/prediction set: a data frame of sample and effect sizes.
#' @param u decision threshold.
#' @param pars list of beta parameter estimates under the alternative (beta1, beta2, sig1)
#' @return data frame of posterior mean and 95% CIs for power in test set
#'
predA = function(testcases, u, pars) {
  psh = mapply(pEstA, theta = testcases$theta, n = testcases$n, MoreArgs = list(u = u, pars = pars))
  df = data.frame(psh = t(psh), n = testcases$n, theta = testcases$theta)
  names(df) = c('psh', 'psl', 'psu', 'n', 'eff')
  df$eff = as.factor(df$eff)
  return(df)
}

#' T1E fit and prediction main function
#'
#' @export
#' @param trcases training set: sample and effect sizes at which simulations are available.
#' @param dat dataframe of simulated sampling distribution in training set.
#' @param testcases test/prediction set: a data frame of sample and effect sizes.
#' @param u decision threshold.
#' @return dataframe of posterior mean and 95% CIs for t1e in test set
#'
T1E = function(trcases, dat) {
  H0index_train = which(trcases$eff==0)
  dfH0 = NULL
  for (j in H0index_train) {
    ps = dat$p_sup[dat$n==trcases$n[j] & dat$eff==trcases$eff[j]]
    a = qmatchH0(ps)
    dfH0 = rbind(dfH0, cbind(n = rep(trcases$n[j], 100), a =a))
  }
  dfH0 = data.frame(dfH0)
  pars = abfith0(dfH0)
  return(pars)
}

#' power fit and prediction main function
#'
#' @export
#' @param trcases training set: sample and effect sizes at which simulations are available.
#' @param dat dataframe of simulated sampling distribution in training set.
#' @param testcases test/prediction set: a data frame of sample and effect sizes.
#' @param u decision threshold.
#' @return dataframe of posterior mean and 95% CIs for power in test set
#'
power = function(trcases, dat) {
  HAindex_train = which(trcases$eff!=0)
  dfHA = NULL
  for (j in HAindex_train) {
    ps = dat$p_sup[dat$n==trcases$n[j] & dat$eff==trcases$eff[j]]
    a = qmatchHA(ps)
    dfHA = rbind(dfHA, cbind(theta = rep(trcases$eff[j], 100), n = rep(trcases$n[j], 100), a =a))
  }
  dfHA = data.frame(dfHA)
  pars = abfitha(dfHA)
  return(pars)
}
