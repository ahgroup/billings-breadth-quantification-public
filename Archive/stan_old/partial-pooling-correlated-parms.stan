//
// Partial Pooling Model with correlated slope and intercept
// Zane
// 2022-02-21
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
// A slope a, and an intercept b. For this model we additionally allow
// each individual to
parameters {
	// Individual-level slope and intercept
	vector[k] a_id;
	vector[k] b_id;

	// Sample-level mean slope and intercept (the partial pooling model will
	// pull outlying lines towards these, which are adaptively estimated).
	real a;
	real b;

	// Variance parameters: s_id are the variances of the individual parameters,
	// s is the variance of the outcome away from the mean, and
	// Rho is the correlation matrix for the slope/intercept.
	vector<lower=0>[2] s_id;
	real<lower=0> s;
	cholesky_factor_corr[2] Lcorr;

	// "Fake" mu which will be drawn from standard normal. Then we can transform
	// it to have the correct distribution -- noncentered parametrization.
	vector[N] mu_tilde;
}

// Backtransform the standard normal mu_tilde to get the mu on the correct
// scale, so we don't have to sample mu directly. Can help avoid "funnel"
// effect observed in multilevel models.
transformed parameters {
	vector[N] mu;
	for (i in 1:N) {
		mu[i] = (a_id[id[i]] + b_id[id[i]] * norm_dist[i]) + s * mu_tilde[i];
	}
}

// The model to be estimated.
// The outcome, y, is normally distributed, where the conditional mean
// on x is fitted to a linear model.
model {
	// Sample parameters that only need to be draw once. These distributions
	// don't depend on anything else.
	Lcorr ~ lkj_corr_cholesky(2);
	s ~ exponential(1);
	s_id ~ exponential(1);
	a ~ normal(0, 1);
	b ~ normal(0, 1);

	// This part deals with the correlated individual-level estimates.
	{
		// Put the population level adaptive prior means into a vector.
		vector[2] mu_pop;
		mu_pop = [a, b]';

		// Vector to hold the individual means. Initialize it for each individual
		// and then draw the values from a multivariate normal distribution.
		array[k] vector[2] mu_id;
		for (j in 1:k) mu_id[j] = [a_id[j], b_id[j]]';
		mu_id ~ multi_normal_cholesky(mu_pop, diag_pre_multiply(s_id, Lcorr));
	}

	// Sample the mu_tilde that will be back transformed
	mu_tilde ~ normal(0, 1);

	// Distribution of the outcome
	y ~ normal(mu, s);
}

// Evaluate the log likelihood of each sample so we can do PSIS or WAIC
generated quantities {
	vector[N] log_lik;
	for (n in 1:N) {
		log_lik[n] = normal_lpdf(y[n] | a_id[id[n]] + b_id[id[n]] * norm_dist[n], s);
	}
	matrix[2, 2] Rho;
	Rho = multiply_lower_tri_self_transpose(Lcorr);
}

// End of program


// Rethinking-generated oxboys example
// data{
	//     int Occasion[234];
	//     vector[234] height;
	//     vector[234] age;
	//     int Subject[234];
	// }
	// parameters{
		//     vector[26] b_sub;
		//     vector[26] a_sub;
		//     real a;
		//     real b;
		//     vector<lower=0>[2] sigma_sub;
		//     real<lower=0> sigma;
		//     corr_matrix[2] Rho;
		// }
		// model{
			//     vector[234] mu;
			//     Rho ~ lkj_corr( 2 );
			//     sigma ~ exponential( 1 );
			//     sigma_sub ~ exponential( 1 );
			//     b ~ normal( 0 , 5 );
			//     a ~ normal( 0 , 5 );
			//     {
				//     vector[2] YY[26];
				//     vector[2] MU;
				// MU = [ a , b ]';
				// for ( j in 1:26 ) YY[j] = [ a_sub[j] , b_sub[j] ]';
				//     YY ~ multi_normal( MU , quad_form_diag(Rho , sigma_sub) );
				//     }
				//     for ( i in 1:234 ) {
					//         mu[i] = a_sub[Subject[i]] + b_sub[Subject[i]] * age[i];
					//     }
					//     height ~ normal( mu , sigma );
					// }
