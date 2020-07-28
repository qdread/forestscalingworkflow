// Fit 2 segment model to richness bins, with FG as random effect
// QDR 27 Jul 2020

data {
	int<lower=0> N;				// Number of bins across all FGs
	int<lower=0> M;				// Number of FGs
	vector<lower=0>[N] x;		// Diameter or light per area
	vector<lower=0>[N] y;		// Richness per bin width
	int<lower=1,upper=M> fg[N];	// Mapping to functional groups 1-M
}

transformed data {
	vector[N] log10_x;
	vector[N] log10_y;
	log10_x = log10(x);
	log10_y = log10(y);
}

parameters {
	// All of the parameters (intercept, slopes, and cutpoints) will have a random effect for FG.
	// Fixed effects:
	real alpha;							// Mean intercept
	real beta_low;						// Mean of low slopes
	real beta_high;						// Mean of high slopes
	real<lower=0> tau;					// Cutpoint
	// Random effects:
	vector[M] alpha_fg;					// Each FG has a deviation from intercept
	vector[M] beta_low_fg;				// Each FG has a deviation from both slopes
	vector[M] beta_high_fg;
	vector[M] tau_fg;					// Each FG has a deviation from the cutpoint
	// Variance parameters
	real<lower=0> sigma;
	real<lower=0,upper=10> sigma_alpha;
	real<lower=0,upper=10> sigma_beta_low;
	real<lower=0,upper=10> sigma_beta_high;
	real<lower=0,upper=10> sigma_tau;
}

transformed parameters {
	vector[N] x2;						// Indicator variable, is x > cutpoint?
	for (i in 1:N) {
		if (log10_x[i] < (tau + tau_fg[fg[i]])) {
			x2[i] = 0;
		} else {
			x2[i] = 1;
		}
	}
}

model {
	// Priors.
	alpha ~ normal(0, 10);
	beta_low ~ normal(0, 10);
	beta_high ~ normal(0, 10);
	tau ~ normal(0, 1);
	alpha_fg ~ normal(0, sigma_alpha);
	beta_low_fg ~ normal(0, sigma_beta_low);
	beta_high_fg ~ normal(0, sigma_beta_high);
	tau_fg ~ normal(0, sigma_tau);
	sigma ~ exponential(1);
	sigma_alpha ~ exponential(1);
	sigma_beta_low ~ exponential(1);
	sigma_beta_high ~ exponential(1);
	sigma_tau ~ exponential(1);
	
	// Likelihood
	log10_y ~ normal((alpha + alpha_fg[fg]) + (beta_low + beta_low_fg[fg]) .* log10_x + (beta_high + beta_high_fg[fg]) .* (log10_x - (tau + tau_fg[fg])) .* x2, sigma);
}

generated quantities {
	// Generate the coefficient estimates for each functional group (fixed+random)
	vector[M] coef_alpha = alpha + alpha_fg;
	vector[M] coef_beta_low = beta_low + beta_low_fg;
	vector[M] coef_beta_high = beta_high + beta_high_fg;
	vector[M] coef_tau = tau + tau_fg;
}

