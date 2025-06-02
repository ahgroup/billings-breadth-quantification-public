//
// Simple linear regression model for simulated data
// WITH CENSORING
// Zane
// 2024-07-18
//

// Input data for the model
// TODO comment this
data {
	// Number of observations
	int<lower=0> N;
	int<lower=0> N_obs;
	int<lower=0> N_cens;

	// Array of observed y values
	array[N_obs] real y_obs;

	// Array of censored y values
	array[N_cens] real y_cens;

	// Predictor vectors -- one where the y's are censored and one where the
	// y's are observed.
	vector[N_obs] x_obs;
	vector[N_cens] x_cens;
}

// The parameters accepted by the model.
// TODO comment this
parameters {
	real alpha;
	real beta;
	real<lower=0> sigma;
}

// TODO comment this
transformed parameters {
	vector[N_obs] mu_obs;
	vector[N_cens] mu_cens;

	mu_obs = alpha + x_obs * beta;
	mu_cens = alpha + x_cens * beta;
}

// The model to be estimated.
// TODO comment this
model {
	// Censored outcome
	y_obs ~ normal(mu_obs, sigma);
	target += normal_lcdf(y_cens | mu_cens, sigma);
}

// End of program
