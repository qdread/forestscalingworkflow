# Extract parameter and fitted values from the mortality fit.

# Function to return fitted value for the logistic regression
logistic_fitted <- function(x, intercept, slope, ...) plogis(intercept + slope * x)

library(tidyverse)
library(rstan)

# Load stan fit
files <- paste0('fit_mortality_', 1:3, '.csv')
fit <- read_stan_csv(file.path('stanoutput', files))

# Calculate summary

##### Get parameters and their credible intervals
sfit <- summary(fit)
params <- data.frame(parameter = row.names(sfit$summary), sfit$summary)

##### Get fitted values and their credible intervals
light_pred <- seq(log10(1), log10(412), length.out = 101)
qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)

# parameters from each sample iteration, then reshape by FG
intercepts <- extract(fit, 'intercept') %>%
	as.data.frame %>%
	setNames(paste0('fg', 1:5)) %>%
	gather(fg, intercept)
slopes <- extract(fit, 'slope') %>%
	as.data.frame %>%
	setNames(paste0('fg', 1:5)) %>%
	gather(fg, slope)

# Combine intercept and slope into one df, then apply function to all
# Get quantiles in same pipe
param_all <- data.frame(iter = 1:3000, intercepts, slope = slopes$slope)	
fitted_quant <- param_all %>%
	group_by(fg, iter) %>%
	group_modify(~ data.frame(light_per_area = 10^light_pred, y = logistic_fitted(x = light_pred, intercept = .x$intercept, slope = .x$slope))) %>%
	group_by(fg, light_per_area) %>%
	group_modify(~ as.data.frame(t(quantile(.x$y, probs = qprobs))) %>% 
		setNames(c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975')))

# Get Bayesian R-squared (Skip for now, can add later)
# To do this the original data are needed for all data points to get the linear predictor.

# Write all output
write_csv(params, 'finalcsvs/mortality_paramci_by_fg.csv')
write_csv(fitted_quant, 'finalcsvs/mortality_ci_by_fg.csv')
