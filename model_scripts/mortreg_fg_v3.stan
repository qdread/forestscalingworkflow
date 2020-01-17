// Logistic regression to model mortality as a function of light received per unit crown area
// Single model fit with random effect for FG (random slope, random intercept model)
// Rewritten to get rid of loops
// QDR 07 Oct 2019
			
data {
	int<lower=0> N;				// Number of trees
	int<lower=0> M;				// Number of FGs
	vector<lower=0>[N] x;		// Light per area
	int<lower=0,upper=1> y[N];	// Mortality 1990-1995
	int<lower=1,upper=M> fg[N];	// Mapping to functional groups 1-M
}

transformed data {
	vector[N] log10_x;
	log10_x = log10(x);
}

parameters {
	real alpha;							// Mean intercept
	vector[M] alpha_fg;					// Deviation of each FG from mean intercept
	real<lower=0,upper=10> sigma_alpha;	// Variation in intercepts
	real beta;							// Mean slope
	vector[M] beta_fg;					// Deviation of each FG from mean slope
	real<lower=0,upper=10> sigma_beta;	// Variation in slopes
	
}

model {
	// Priors
	alpha ~ normal(0, 10);
	beta ~ normal(0, 10);
	alpha_fg ~ normal(0, sigma_alpha);
	beta_fg ~ normal(0, sigma_beta);
	sigma_alpha ~ exponential(1);
	sigma_beta ~ exponential(1);
		
	// Likelihood
	y ~ bernoulli_logit(alpha + alpha_fg[fg] + (beta + beta_fg[fg]) .* log10_x);

}

generated quantities {
	// Generate the coefficient estimates for each functional group (fixed+random)
	vector[M] intercept = alpha + alpha_fg;
	vector[M] slope = beta + beta_fg;
}
