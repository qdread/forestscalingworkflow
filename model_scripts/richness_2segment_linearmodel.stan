// Fit 2 segment model to richness bins, to be done for all bins together (no random effect)
// QDR 27 Jul 2020

data {

	int<lower=0> N;				// Number of bins
	vector<lower=0>[N] x;		// Diameter or light per area
	vector<lower=0>[N] y;		// Richness per bin width
}

transformed data {
	vector[N] log10_x;
	vector[N] log10_y;
	log10_x = log10(x);
	log10_y = log10(y);
}

parameters {
	real alpha;							// Mean intercept
	vector[2] beta;						// Beta low and beta high
	real<lower=0> sigma;				// Variation in intercepts
	real<lower=0> tau;					// Cutpoint
}

transformed parameters {
	vector[N] x2;						// Indicator variable, is x > cutpoint?
	for (i in 1:N) {
		if (log10_x[i] < tau) {
			x2[i] = 0;
		} else {
			x2[i] = 1;
		}
	}
}

model {
	// Priors.
	vector[N] mu;
	alpha ~ normal(0, 1);
	beta ~ normal(0, 1);
	sigma ~ normal(0, 2);
	tau ~ normal(0, 1);
	
	// Likelihood
	for (i in 1:N) {
		mu[i] = alpha + beta[1] * log10_x[i] + beta[2] * (log10_x[i] - tau) * x2[i];
	}
	
	log10_y ~ normal(mu, sigma);
}
