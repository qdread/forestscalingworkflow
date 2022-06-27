# Get credible intervals for ratio slopes, and fitted values for the ratio slope trend (for plotting); for diameter only.
# QDR / Forestlight / 26 Feb 2020

get_ratio_slopes <- function(dens_model, prod_model, fg_top, fg_bottom, year, xmin, n, use_subset = FALSE, n_chains = 3, scaling_var = 'dbh', fp = '~/forestlight/stanoutput', densityfitprefix = 'fit_density', productionfitprefix = 'fit_production', scalingtype = 'production') {
  require(rstan)
  require(Brobdingnag)
  require(purrr)
  
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) 
  # Uses formula x/y dy/dx for log slope.
  log_slope <- function(x, y) (x[-1]/y[-1]) * diff(y)/diff(x)  
  
  # names of parameters
  density_par <- list('1' = c('alpha'),
                      '2' = c('alpha_low', 'alpha_high', 'tau'),
                      '3' = c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high'))
  production_par <- list('1' = c('beta0', 'beta1', 'sigma'),
                         '2' = c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta', 'sigma'))					  
  
  pars_to_get <- list(density = density_par[[dens_model]], production = production_par[[prod_model]])
  
  # Load CSVs as stanfit object
  message('Loading stan fit . . .')
  files_density_top <- paste0(densityfitprefix, dens_model, '_', 'production', '_', fg_top, '_', year, '_', 1:n_chains, '.csv')
  files_production_top <- paste0(productionfitprefix, prod_model, '_', scalingtype, '_', fg_top, '_', year, '_', 1:n_chains, '.csv')
  files_density_bottom <- paste0(densityfitprefix, dens_model, '_', 'production', '_', fg_bottom, '_', year, '_', 1:n_chains, '.csv')
  files_production_bottom <- paste0(productionfitprefix, prod_model, '_', scalingtype, '_', fg_bottom, '_', year, '_', 1:n_chains, '.csv')
  fit <- list(density_top = read_stan_csv(file.path(fp,files_density_top)), 
              production_top = read_stan_csv(file.path(fp,files_production_top)),
              density_bottom = read_stan_csv(file.path(fp,files_density_bottom)), 
              production_bottom = read_stan_csv(file.path(fp,files_production_bottom)))
  
  message('Calculating density slopes . . .')
  # extract density parameters
  pars_dens_top <- as.data.frame(do.call('cbind', extract(fit[['density_top']], pars_to_get[['density']])))
  pars_dens_bottom <- as.data.frame(do.call('cbind', extract(fit[['density_bottom']], pars_to_get[['density']])))
  
  # get fitted density values
  if (dens_model == '1') {
    dens_fitted_top <- sapply(dbh_pred, pdf_pareto, xmin = xmin, alpha = pars_dens_top[,'alpha']) 
    dens_fitted_bottom <- sapply(dbh_pred, pdf_pareto, xmin = xmin, alpha = pars_dens_bottom[,'alpha']) 
  }
  if (dens_model == '2') {
    dens_fitted_top <- sapply(dbh_pred, pdf_2part, xmin = xmin, alpha_low = pars_dens_top[,'alpha_low'], alpha_high = pars_dens_top[,'alpha_high'], tau = pars_dens_top[,'tau'])
    dens_fitted_bottom <- sapply(dbh_pred, pdf_2part, xmin = xmin, alpha_low = pars_dens_bottom[,'alpha_low'], alpha_high = pars_dens_bottom[,'alpha_high'], tau = pars_dens_bottom[,'tau'])
  }
  if (dens_model == '3') {
    dens_fitted_top <- sapply(dbh_pred, pdf_3part, xmin = xmin, alpha_low = pars_dens_top[,'alpha_low'], alpha_mid = pars_dens_top[,'alpha_mid'], alpha_high = pars_dens_top[,'alpha_high'], tau_low = pars_dens_top[,'tau_low'], tau_high = pars_dens_top[,'tau_high'])
    dens_fitted_bottom <- sapply(dbh_pred, pdf_3part, xmin = xmin, alpha_low = pars_dens_bottom[,'alpha_low'], alpha_mid = pars_dens_bottom[,'alpha_mid'], alpha_high = pars_dens_bottom[,'alpha_high'], tau_low = pars_dens_bottom[,'tau_low'], tau_high = pars_dens_bottom[,'tau_high'])
  }
  
  # calculate fitted slopes
  dens_fitted_top <- dens_fitted_top * n[1]
  dens_fittedslopes_top <- map_dfr(as.data.frame(t(dens_fitted_top)), ~ log_slope(dbh_pred, .))
  
  dens_fitted_bottom <- dens_fitted_bottom * n[2]
  dens_fittedslopes_bottom <- map_dfr(as.data.frame(t(dens_fitted_bottom)), ~ log_slope(dbh_pred, .))
  
  message('Getting production values . . .')
  # extract production parameters
  pars_prod_top <- as.data.frame(do.call('cbind', extract(fit[['production_top']], pars_to_get[['production']])))
  pars_prod_bottom <- as.data.frame(do.call('cbind', extract(fit[['production_bottom']], pars_to_get[['production']])))
  
  # get fitted production values
  if (prod_model == '1') {
    prod_fitted_top <- sapply(dbh_pred, powerlaw_log, beta0 = pars_prod_top[,'beta0'], beta1 = pars_prod_top[,'beta1'])
    prod_fitted_bottom <- sapply(dbh_pred, powerlaw_log, beta0 = pars_prod_bottom[,'beta0'], beta1 = pars_prod_bottom[,'beta1'])
  }
  if (prod_model == '2') {
    prod_fitted_top <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars_prod_top[,'beta0'], beta1_low = pars_prod_top[,'beta1_low'], beta1_high = pars_prod_top[,'beta1_high'], x0 = pars_prod_top[,'x0'], delta = pars_prod_top[,'delta'])
    prod_fitted_bottom <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars_prod_bottom[,'beta0'], beta1_low = pars_prod_bottom[,'beta1_low'], beta1_high = pars_prod_bottom[,'beta1_high'], x0 = pars_prod_bottom[,'x0'], delta = pars_prod_bottom[,'delta'])
  }
  
  message('Calculating total production slopes  . . .')
  
  totalprod_fitted_top <- dens_fitted_top * prod_fitted_top
  totalprod_fittedslopes_top <- map_dfr(as.data.frame(t(totalprod_fitted_top)), ~ log_slope(dbh_pred, .))
  
  totalprod_fitted_bottom <- dens_fitted_bottom * prod_fitted_bottom
  totalprod_fittedslopes_bottom <- map_dfr(as.data.frame(t(totalprod_fitted_bottom)), ~ log_slope(dbh_pred, .))
  
  message('Getting quantiles of ratio slopes (almost done!) . . .')
  # Ratio of density slopes
  ratio_fittedslopes_dens <- dens_fittedslopes_top - dens_fittedslopes_bottom
  # Ratio of total production slopes
  ratio_fittedslopes_totalprod <- totalprod_fittedslopes_top - totalprod_fittedslopes_bottom
  
  # Quantiles for each ratio
  ratio_fittedslopes_dens_quant <- apply(ratio_fittedslopes_dens, 1, quantile, probs = qprobs, na.rm = TRUE)
  ratio_fittedslopes_totalprod_quant <- apply(ratio_fittedslopes_totalprod, 1, quantile, probs = qprobs, na.rm = TRUE)
  
  # Process output into neat data frame and return
  out <- data.frame(dbh = dbh_pred[-1],
                    variable = rep(c('density', 'total production'), each = length(dbh_pred) - 1),
                    rbind(t(ratio_fittedslopes_dens_quant), t(ratio_fittedslopes_totalprod_quant)))
  
  setNames(out, c(scaling_var, 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}

# Fitted value function also takes input of the production data for correction factor calculations.
get_ratio_fitted <- function(dens_model, prod_model, prod_data, fg_top, fg_bottom, year, xmin, n, use_subset = FALSE, n_chains = 3, scaling_var = 'dbh', fp = '~/forestlight/stanoutput', densityfitprefix = 'fit_density', productionfitprefix = 'fit_production', scalingtype = 'production') {
  require(rstan)
  require(Brobdingnag)
  require(purrr)
  
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) 
  # Uses formula x/y dy/dx for log slope.
  log_slope <- function(x, y) (x[-1]/y[-1]) * diff(y)/diff(x)  
  
  # names of parameters
  density_par <- list('1' = c('alpha'),
                      '2' = c('alpha_low', 'alpha_high', 'tau'),
                      '3' = c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high'))
  production_par <- list('1' = c('beta0', 'beta1', 'sigma'),
                         '2' = c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta', 'sigma'))					  
  
  pars_to_get <- list(density = density_par[[dens_model]], production = production_par[[prod_model]])
  
  # Load CSVs as stanfit object
  message('Loading stan fit . . .')
  files_density_top <- paste0(densityfitprefix, dens_model, '_', 'production', '_', fg_top, '_', year, '_', 1:n_chains, '.csv')
  files_production_top <- paste0(productionfitprefix, prod_model, '_', scalingtype, '_', fg_top, '_', year, '_', 1:n_chains, '.csv')
  files_density_bottom <- paste0(densityfitprefix, dens_model, '_', 'production', '_', fg_bottom, '_', year, '_', 1:n_chains, '.csv')
  files_production_bottom <- paste0(productionfitprefix, prod_model, '_', scalingtype, '_', fg_bottom, '_', year, '_', 1:n_chains, '.csv')
  fit <- list(density_top = read_stan_csv(file.path(fp,files_density_top)), 
              production_top = read_stan_csv(file.path(fp,files_production_top)),
              density_bottom = read_stan_csv(file.path(fp,files_density_bottom)), 
              production_bottom = read_stan_csv(file.path(fp,files_production_bottom)))
  
  message('Getting density fitted values . . .')
  # extract density parameters
  pars_dens_top <- as.data.frame(do.call('cbind', extract(fit[['density_top']], pars_to_get[['density']])))
  pars_dens_bottom <- as.data.frame(do.call('cbind', extract(fit[['density_bottom']], pars_to_get[['density']])))
  
  # get fitted density values
  if (dens_model == '1') {
    dens_fitted_top <- sapply(dbh_pred, pdf_pareto, xmin = xmin, alpha = pars_dens_top[,'alpha']) 
    dens_fitted_bottom <- sapply(dbh_pred, pdf_pareto, xmin = xmin, alpha = pars_dens_bottom[,'alpha']) 
  }
  if (dens_model == '2') {
    dens_fitted_top <- sapply(dbh_pred, pdf_2part, xmin = xmin, alpha_low = pars_dens_top[,'alpha_low'], alpha_high = pars_dens_top[,'alpha_high'], tau = pars_dens_top[,'tau'])
    dens_fitted_bottom <- sapply(dbh_pred, pdf_2part, xmin = xmin, alpha_low = pars_dens_bottom[,'alpha_low'], alpha_high = pars_dens_bottom[,'alpha_high'], tau = pars_dens_bottom[,'tau'])
  }
  if (dens_model == '3') {
    dens_fitted_top <- sapply(dbh_pred, pdf_3part, xmin = xmin, alpha_low = pars_dens_top[,'alpha_low'], alpha_mid = pars_dens_top[,'alpha_mid'], alpha_high = pars_dens_top[,'alpha_high'], tau_low = pars_dens_top[,'tau_low'], tau_high = pars_dens_top[,'tau_high'])
    dens_fitted_bottom <- sapply(dbh_pred, pdf_3part, xmin = xmin, alpha_low = pars_dens_bottom[,'alpha_low'], alpha_mid = pars_dens_bottom[,'alpha_mid'], alpha_high = pars_dens_bottom[,'alpha_high'], tau_low = pars_dens_bottom[,'tau_low'], tau_high = pars_dens_bottom[,'tau_high'])
  }
  
  # calculate fitted values
  dens_fitted_top <- dens_fitted_top * n[1]
  
  dens_fitted_bottom <- dens_fitted_bottom * n[2]
  
  message('Getting production fitted values . . .')
  # extract production parameters
  pars_prod_top <- as.data.frame(do.call('cbind', extract(fit[['production_top']], pars_to_get[['production']])))
  pars_prod_bottom <- as.data.frame(do.call('cbind', extract(fit[['production_bottom']], pars_to_get[['production']])))
  
  # get fitted production values
  if (prod_model == '1') {
    prod_fitted_top <- sapply(dbh_pred, powerlaw_log, beta0 = pars_prod_top[,'beta0'], beta1 = pars_prod_top[,'beta1'])
    prod_fitted_bottom <- sapply(dbh_pred, powerlaw_log, beta0 = pars_prod_bottom[,'beta0'], beta1 = pars_prod_bottom[,'beta1'])
  }
  if (prod_model == '2') {
    prod_fitted_top <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars_prod_top[,'beta0'], beta1_low = pars_prod_top[,'beta1_low'], beta1_high = pars_prod_top[,'beta1_high'], x0 = pars_prod_top[,'x0'], delta = pars_prod_top[,'delta'])
    prod_fitted_bottom <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars_prod_bottom[,'beta0'], beta1_low = pars_prod_bottom[,'beta1_low'], beta1_high = pars_prod_bottom[,'beta1_high'], x0 = pars_prod_bottom[,'x0'], delta = pars_prod_bottom[,'delta'])
  }
  
  message('Calculating total production correction factor . . .')
  # Get total production including correction factors.

	corr_factor <- function(y, y_fit, n_pars) {
	  y_fit <- do.call(cbind, y_fit)
	  # Get residuals by subtracting log y from linear predictor
	  resids <- -1 * sweep(log(y_fit), 1, log(y))
	  
	  # Sum of squared residuals
	  ssq_resid <- apply(resids^2, 2, sum)
	  # Standard error of estimates
	  sse <- (ssq_resid / (length(y) - n_pars))^0.5
	  # Correction factors
	  exp((sse^2)/2)
	}

	# We need all fitted values for production, not just the 101 values, to calculate correction factor.
	prod_pred_all_top <- pmap(pars_prod_top[, c('beta0','beta1')], powerlaw_log, x = prod_data[[1]]$x)
	prod_pred_all_bottom <- pmap(pars_prod_bottom[, c('beta0','beta1')], powerlaw_log, x = prod_data[[2]]$x)
	
	prod_cf_top <- corr_factor(y = prod_data[[1]]$y, y_fit = prod_pred_all_top, n_pars = 2)
	prod_cf_bottom <- corr_factor(y = prod_data[[2]]$y, y_fit = prod_pred_all_bottom, n_pars = 2)

  
  
  # get fitted total production values
  totalprod_fitted_top <- sweep(dens_fitted_top * prod_fitted_top, 1, prod_cf_top, `*`)
  totalprod_fitted_bottom <- sweep(dens_fitted_bottom * prod_fitted_bottom, 1, prod_cf_bottom, `*`)
  
  message('Getting quantiles of fitted values (almost done!) . . .')
  # Ratio of density fitted values
  ratio_dens <- dens_fitted_top / dens_fitted_bottom
  # Ratio of total production fitted values
  ratio_totalprod <- totalprod_fitted_top / totalprod_fitted_bottom
  
  # Quantiles for each ratio
  ratio_dens_quant <- apply(ratio_dens, 2, quantile, probs = qprobs, na.rm = TRUE)
  ratio_totalprod_quant <- apply(ratio_totalprod, 2, quantile, probs = qprobs, na.rm = TRUE)
  
  # Process output into neat data frame and return
  out <- data.frame(dbh = dbh_pred,
                    variable = rep(c('density', 'total production'), each = length(dbh_pred)),
                    rbind(t(ratio_dens_quant), t(ratio_totalprod_quant)))
  
  setNames(out, c(scaling_var, 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}


# Execute function for all pairwise ratios --------------------------------

# Modified 27 June 2022: this is now done for all six pairwise combinations of the four named FGs instead of just selected ones.

dbh_pred <- exp(seq(log(1), log(315), length.out = 101))
fg_pairs <- combn(1:4, 2, simplify = FALSE)
fg_titles <- c('fast', 'tall', 'slow', 'short')

library(forestscaling)
library(dplyr)
library(purrr)

min_n <- read.csv('~/forestlight/stanrdump/min_n.csv')

# Get slopes and credible intervals
ratioslopes <- map_dfr(fg_pairs, function(fgs) {
  fg_names <- paste0('fg', fgs)
  ns <- min_n$n[min_n$fg %in% fg_names & min_n$year == 1995]
  slopes <- get_ratio_slopes(dens_model = '3', prod_model = '1', fg_top = fg_names[1], fg_bottom = fg_names[2], year = 1995, xmin = 1, n = ns)
  data.frame(ratio = paste(fg_titles[fgs], collapse = ':'),
             slopes)
})

write.csv(ratioslopes, '~/forestlight/finalcsvs/ratio_slope_ci_allpairwise.csv', row.names = FALSE)

# Fitted values

# load data dumps.
prod_dumps <- map(1:4, function(i) {
	dumpfile <- paste0('dump_production_fg', i, '_1995.r')
	source(file.path('~/forestlight/stanrdump', dumpfile)) # Creates variables x and y.
	list(x = x, y = y)
})

ratiofitted <- map_dfr(fg_pairs, function(fgs) {
  fg_names <- paste0('fg', fgs)
  ns <- min_n$n[min_n$fg %in% fg_names & min_n$year == 1995]
  fittedvals <- get_ratio_fitted(dens_model = '3', prod_model = '1', prod_data = prod_dumps[fgs], fg_top = fg_names[1], fg_bottom = fg_names[2], year = 1995, xmin = 1, n = ns)
  data.frame(ratio = paste(fg_titles[fgs], collapse = ':'),
             fittedvals)
})

write.csv(ratiofitted, '~/forestlight/finalcsvs/ratio_fittedvalues_allpairwise.csv', row.names = FALSE)
