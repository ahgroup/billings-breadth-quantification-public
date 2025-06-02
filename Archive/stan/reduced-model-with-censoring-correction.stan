// generated with brms 2.21.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  array[N] int<lower=-1,upper=2> cens;  // indicates censoring
  vector[N] rcens;  // right censor points for interval censoring
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept | 3, 0, 3);
  lprior += student_t_lpdf(sigma | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;
    for (n in 1:N) {
    // special treatment of censored data
      if (cens[n] == 0) {
        target += normal_lpdf(Y[n] | mu[n], sigma);
      } else if (cens[n] == 1) {
        target += normal_lccdf(Y[n] | mu[n], sigma);
      } else if (cens[n] == -1) {
        target += normal_lcdf(Y[n] | mu[n], sigma);
      } else if (cens[n] == 2) {
        target += log_diff_exp(
          normal_lcdf(rcens[n] | mu[n], sigma),
          normal_lcdf(Y[n] | mu[n], sigma)
        );
      }
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}

