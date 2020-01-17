// Bertalanffy model with log link function
// Used to model production per unit area for individual trees as a function of incoming light per unit area
// New version created on 15 June with correct way of getting the max slope
			
data {
	int<lower=0> N;
	vector<lower=0>[N] x;	// Light per area
	vector<lower=0>[N] y;	// Production per area
}

transformed data {
	vector[N] log_y;
	log_y = log(y);
}

parameters {
	real<lower=0> A;
	real<lower=0> b;
	real<lower=0> k;
	real<lower=0> sigma;
}

model {
	// Priors
	A ~ lognormal(1, 1);
	b ~ lognormal(1, 1);
	k ~ lognormal(1, 1);
		
	// Likelihood
	{
		vector[N] mu;
		for (i in 1:N) mu[i] = A * (1 - b * exp(-k * x[i])) ^ 3;
		log_y ~ normal(log(mu), sigma);
	}
}

generated quantities {
	real x_max;
	real y_max;
	real log_slope;
	
	// Analytical formulas added 15 Jun 2018
	x_max = (1 - b/exp(1.0)) / k;
	y_max = A * (1 - b * exp(-k * x_max)) ^ 3;
	log_slope = 3*b*k * (x_max*exp(-k*x_max)) / (1 - b*exp(-k*x_max));
} 
