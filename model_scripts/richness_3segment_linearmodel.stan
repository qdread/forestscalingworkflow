// Fit 3 segment model to richness bins, to be done for all bins together (no random effect)
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
	real alpha;							// Intercept
	vector[3] beta;						// Beta low and beta high
	real<lower=0> sigma;				// Variance
	real<lower=0, upper=1> tau_low;		// Cutpoint 1. Constrain to be low
	real<lower=1.3, upper=3> tau_high;	// Cutpoint 2. Constrain to be sufficiently greater than cutpoint 1.
}

transformed parameters {
	vector[N] x2a;						// Indicator variable, is x > cutpoint 1 but <= cutpoint 2
	vector[N] x2b;						// Indicator variable, is x > cutpoint 2
	for (i in 1:N) {
		if (log10_x[i] < tau_low) {
			x2a[i] = 0;
		} else {
			x2a[i] = 1;
		}
		if (log10_x[i] < tau_high) {
			x2b[i] = 0;
		} else {
			x2b[i] = 1;
		}
	}
}

model {
	// Priors.
	vector[N] mu;
	alpha ~ normal(0, 1);
	beta ~ normal(0, 1);
	sigma ~ normal(0, 2);
	// Taus have no prior (uniform on supported interval)
	
	// Likelihood
	for (i in 1:N) {
		mu[i] = alpha + beta[1] * log10_x[i] + beta[2] * (log10_x[i] - tau_low) * x2a[i] + beta[3] * (log10_x[i] - tau_high) * x2b[i];
	}
	
	log10_y ~ normal(mu, sigma);
}
