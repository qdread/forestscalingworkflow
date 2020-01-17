data {
	int<lower=0> N;
	vector<lower=0>[N] x;
	vector<lower=0>[N] y;
    real<lower=0> x_min;
    real<lower=0> x_max;
}

transformed data {
	vector[N] logx;
	logx = log(x);
}

parameters {
	// Pareto density
	real<lower=0, upper=5> alpha;
}

model {
	// Prior: Pareto density
	alpha ~ lognormal(1, 1) T[0, 5];
	
	// Likelihood: Pareto density
	x ~ pareto(x_min, alpha);
}

generated quantities {
	vector[N] log_lik; // Log-likelihood for getting info criteria later
	
	for (i in 1:N) {
		log_lik[i] = pareto_lpdf(x[i] | x_min, alpha);
	}
}
