//
// No Pooling Model where season is ignored
// Zane
// 2022-02-18
//

// Input data for the model
// Input data is of length N.
// There are two variables that we will use: y and x.
data {
  int<lower=0> N;
  int<lower=0> k;
  array[N] int id;
  vector[N] y;
  vector[N] norm_dist;
}

// The parameters accepted by the model.
// Our model has a variance parameter, s, and two linear model parameters.
// A slope a, and an intercept b.
parameters {
  vector[k] a;
  vector[k] b;
  real<lower=0> s;
}

// The model to be estimated.
// The outcome, y, is normally distributed, where the conditional mean
// on x is fitted to a linear model.
model {
	vector[N] mu;
	s ~ exponential(1);
	b ~ normal(0, 2);
	a ~ normal(5, 2);
	for (i in 1:N) {
		mu[i] = a[id[i]] + b[id[i]] * norm_dist[i];
	};
  y ~ normal(mu, s);
}

// Evaluate the log likelihood of each sample so we can do PSIS or WAIC
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] |  a[id[n]] + b[id[n]] * norm_dist[n], s);
  }

  // Prior predictive check
  // vector[N] ones_N = rep_vector(1, N);
  // vector[N] a_p = normal_rng(ones_N * 5, ones_N * 2);
  // vector[N] b_p = normal_rng(ones_N * 0, ones_N * 2);
  // real s_p = exponential_rng(1);
  // array[N] real y_prior = normal_rng(a_p + b_p * x, s_p);
}

// End of program
