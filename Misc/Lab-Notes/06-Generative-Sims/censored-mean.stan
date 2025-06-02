//
// Estimating mean of censored data NO DSITANCES!!!!
// Zane
// 2024-07-18
//

// Input data for the model
// Input data is of length N.
// There are two variables that we will use: y and x.
data {
	// Different count variables
	int<lower=0> N_obs;
	int<lower=0> N_cens;

	// Array of observed y values
	array[N_obs] real y_obs;

	// Point below which data is left-censored
	real<upper=min(y_obs)> U;
}

// The parameters accepted by the model.
parameters {
	array[N_cens] real<upper=U> y_cens;
	real mu;
	real<lower=0> sigma;
}

// The model to be estimated.
// The outcome, y, is normally distributed, where the conditional mean
// on x is fitted to a linear model.
model {
	// Censored outcome
	y_obs ~ normal(mu, sigma);
	y_cens ~ normal(mu, sigma);
}

// End of program
