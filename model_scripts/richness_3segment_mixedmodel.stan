// Fit 3 segment model to richness bins, with FG as random effect
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
	real beta_mid;						// Mean of mid slopes
	real beta_high;						// Mean of high slopes
	real<lower=0> tau_low;				// Cutpoint
	real<lower=0> tau_distance;			// Distance between tau low and tau high, must > 0
	// Random effects:
	vector[M] alpha_fg;					// Each FG has a deviation from intercept
	vector[M] beta_low_fg;				// Each FG has a deviation from all slopes
	vector[M] beta_mid_fg;
	vector[M] beta_high_fg;
	vector[M] tau_low_fg;				// Each FG has a deviation from the cutpoint
	vector<lower=-1*tau_distance>[M] tau_distance_fg;
	// Variance parameters
	real<lower=0> sigma;
	real<lower=0,upper=10> sigma_alpha;
	real<lower=0,upper=10> sigma_beta_low;
	real<lower=0,upper=10> sigma_beta_mid;
	real<lower=0,upper=10> sigma_beta_high;
	real<lower=0,upper=10> sigma_tau_low;
	real<lower=0,upper=10> sigma_tau_distance;
}

transformed parameters {
	real<lower=0> tau_high;
	vector[M] tau_high_fg;
	vector[N] x2a;						// Indicator variable, is x > cutpoint?
	vector[N] x2b;
	
	tau_high = tau_low + tau_distance;
	tau_high_fg = tau_low_fg + tau_distance_fg;
	
	for (i in 1:N) {
		if (log10_x[i] < (tau_low + tau_low_fg[fg[i]])) {
			x2a[i] = 0;
		} else {
			x2a[i] = 1;
		}
		if (log10_x[i] < (tau_high + tau_high_fg[fg[i]])) {
			x2b[i] = 0;
		} else {
			x2b[i] = 1;
		}
	}
}

model {
	// Priors.
	alpha ~ normal(0, 10);
	beta_low ~ normal(0, 10);
	beta_mid ~ normal(0, 10);
	beta_high ~ normal(0, 10);
	tau_low ~ normal(0, 1);
	tau_distance ~ normal(0, 1);
	alpha_fg ~ normal(0, sigma_alpha);
	beta_low_fg ~ normal(0, sigma_beta_low);
	beta_mid_fg ~ normal(0, sigma_beta_mid);
	beta_high_fg ~ normal(0, sigma_beta_high);
	tau_low_fg ~ normal(0, sigma_tau_low);
	tau_distance_fg ~ normal(0, sigma_tau_distance);
	sigma ~ exponential(1);
	sigma_alpha ~ exponential(1);
	sigma_beta_low ~ exponential(1);
	sigma_beta_mid ~ exponential(1);
	sigma_beta_high ~ exponential(1);
	sigma_tau_low ~ exponential(1);
	sigma_tau_distance ~ exponential(1);
	
	// Likelihood
	log10_y ~ normal((alpha + alpha_fg[fg]) + (beta_low + beta_low_fg[fg]) .* log10_x + (beta_mid + beta_mid_fg[fg]) .* (log10_x - (tau_low + tau_low_fg[fg])) .* x2a + (beta_high + beta_high_fg[fg]) .* (log10_x - (tau_high + tau_high_fg[fg])) .* x2b, sigma);
}

generated quantities {
	// Generate the coefficient estimates for each functional group (fixed+random)
	vector[M] coef_alpha = alpha + alpha_fg;
	vector[M] coef_beta_low = beta_low + beta_low_fg;
	vector[M] coef_beta_mid = beta_mid + beta_mid_fg;
	vector[M] coef_beta_high = beta_high + beta_high_fg;
	vector[M] coef_tau_low = tau_low + tau_low_fg;
	vector[M] coef_tau_high = tau_low + tau_low_fg + tau_distance + tau_distance_fg;
}

