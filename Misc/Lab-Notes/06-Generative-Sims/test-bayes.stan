//
// Simple linear regression model for simulated data
// Zane
// 2023-06-13
//

// Input data for the model
// Input data is of length N.
// There are two variables that we will use: y and x.
data {
  int<lower=0> N;
  real<lower=0> s_mean;
  real a_mean;
  real b_mean;
  real<lower=0> a_sd;
  real<lower=0> b_sd;
  vector[N] y;
  vector[N] d;
}

// The parameters accepted by the model.
// Our model has a variance parameter, s, and two linear model parameters.
// A slope a, and an intercept b.
parameters {
  real a;
  real b;
  real<lower=0> s;
}

// The model to be estimated.
// The outcome, y, is normally distributed, where the conditional mean
// on x is fitted to a linear model.
model {
	vector[N] mu;
	s ~ exponential(s_mean);
	b ~ normal(b_mean, b_sd);
	a ~ normal(a_mean, a_sd);
	for (i in 1:N) {
		mu[i] = a + b * d[i];
	}
  y ~ normal(mu, s);
}

// End of program
