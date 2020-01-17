# Source code for all functions to extract information from piecewise model fits
# Part of the final stan workflow
# QDR / Forestlight / 10 Sep 2018
# Edited 14 Feb 2019: also expand this to handle the piecewise by light received fits.
# Edited 17 Feb 2019: Add lognormal distribution too.
# Edited 17 Jun 2019: Change file names so it can accommodate the incoming light and volume scalings
# Edited 02 Jul 2019: Get rid of the integration correction factor (better new one is applied elsewhere)
# Edited 10 Sep 2019: Also extract the sigma parameter from production fits.

# Information to get:
# Parameter values and credible intervals
# Fitted values and credible intervals (and prediction intervals)
# Fitted log slopes and credible intervals
# Information criteria
# "Bayesian" R-squared
# Corrected hinge function. Must accommodate very high numbers.

# 0. Define all density and production functions.
# -----------------------------------------------

# Now in package.
library(forestscaling)

# 1. Function for parameter values and their credible intervals
# -------------------------------------------------------------

param_values <- function(fit, pars) {
   summ_fit <- summary(fit)
   param_cis <- summ_fit[[1]][pars, , drop = FALSE]

   param_cis <- data.frame(parameter = dimnames(param_cis)[[1]], param_cis)
   names(param_cis)[5:9] <- c('q025', 'q25', 'q50', 'q75', 'q975')

   return(param_cis)
}

# 2. Functions to get fitted values and credible/prediction intervals
# ------------------------------------------------------------------

# Works for density, for production, or for total production
# Fit should be a list. If total production, it's list (density, production)

fitted_predicted_values <- function(fit, variable, dbh_pred, dens_form = NA, prod_form = NA, x_min = NULL, n_indiv = 1, ll = 1.1, ul = 316, pars_to_get, scaling_var = 'dbh') {
  require(purrr)

  # Use high and low normal quantiles to get prediction intervals
  prediction_quantile <- function(p, fitted_values) {
    pred_raw <- apply(fitted_values, 2, function(mus) exp(qnorm(p, mean = log(mus), sd = pars_prod[,'sigma'])))
    apply(pred_raw, 2, median, na.rm = TRUE)
  }
  
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) 
  
  if (variable == 'density' | variable == 'total_production') {
    pars_dens <- as.data.frame(do.call('cbind', extract(fit[['density']], pars_to_get[['density']])))
    
    if (dens_form == '1') {
      dens_fitted <- sapply(dbh_pred, pdf_pareto, xmin = x_min, alpha = pars_dens[,'alpha']) 
    }
    if (dens_form == '2') {
      dens_fitted <- sapply(dbh_pred, pdf_2part, xmin = x_min, alpha_low = pars_dens[,'alpha_low'], alpha_high = pars_dens[,'alpha_high'], tau = pars_dens[,'tau'])
    }
    if (dens_form == '3') {
      dens_fitted <- sapply(dbh_pred, pdf_3part, xmin = x_min, alpha_low = pars_dens[,'alpha_low'], alpha_mid = pars_dens[,'alpha_mid'], alpha_high = pars_dens[,'alpha_high'], tau_low = pars_dens[,'tau_low'], tau_high = pars_dens[,'tau_high'])
    }
    if (dens_form == 'w') {
      # Must manually rescale and remove upper and lower truncations
      trunc_pts <- pmap(pars_dens, function(m, n, ...) pweibull(q = c(ll,ul), shape = m, scale = n))
      dens_fitted <- sapply(dbh_pred, dweibull, shape = pars_dens[,'m'], scale = pars_dens[,'n']) 
      dens_fitted <- t(sapply(1:nrow(dens_fitted), function(i) {
        x <- dens_fitted[i,]/diff(trunc_pts[[i]])
        x
      }))
    }
    if (dens_form == 'ln') {
      dens_fitted <- sapply(dbh_pred, dlnorm, meanlog = pars_dens[,'mu_logn'], sdlog = pars_dens[,'sigma_logn']) 
    }
    dens_fitted <- dens_fitted * n_indiv
    dens_fitted_quant <- apply(dens_fitted, 2, quantile, probs = qprobs, na.rm = TRUE)
  }
  
  if (variable == 'production' | variable == 'total_production') {
    pars_prod <- as.data.frame(do.call('cbind', extract(fit[['production']], c(pars_to_get[['production']]))))
    if (prod_form == '1') {
      prod_fitted <- sapply(dbh_pred, powerlaw_log, beta0 = pars_prod[,'beta0'], beta1 = pars_prod[,'beta1'])
    }
    if (prod_form == '2') {
      prod_fitted <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars_prod[,'beta0'], beta1_low = pars_prod[,'beta1_low'], beta1_high = pars_prod[,'beta1_high'], x0 = pars_prod[,'x0'], delta = pars_prod[,'delta'])
    }
    prod_fitted_quant <- apply(prod_fitted, 2, quantile, probs = qprobs, na.rm = TRUE)
    prod_pred_quant <- map_dfc(qprobs, prediction_quantile, fitted_values = prod_fitted)
  }
  
  if (variable == 'total_production') {
    totalprod_fitted <- dens_fitted * prod_fitted
    
    totalprod_fitted_quant <- apply(totalprod_fitted, 2, quantile, probs = qprobs, na.rm = TRUE)
    totalprod_pred_quant <- map_dfc(qprobs, prediction_quantile, fitted_values = totalprod_fitted)
  }
  
  if (variable == 'density') {
    out <- data.frame(dbh = dbh_pred,
                      variable = rep('density', length(dbh_pred)),
                      t(dens_fitted_quant))
  } 
  if (variable == 'production') {
    out <- data.frame(dbh = dbh_pred,
                      variable = rep(c('production','production_fitted'), each = length(dbh_pred)),
                      rbind(as.matrix(prod_pred_quant), t(prod_fitted_quant)))
  }
  if (variable == 'total_production') {
    out <- data.frame(dbh = dbh_pred,
                      variable = rep(c('total_production','total_production_fitted'), each = length(dbh_pred)),
                      rbind(as.matrix(totalprod_pred_quant), t(totalprod_fitted_quant)))
  }
  
  setNames(out, c(scaling_var, 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))  
}


# 3. Function to get the fitted log slope values at many points
# -------------------------------------------------------------

fitted_slope_ci <- function(fit, variable, dbh_pred, dens_form, prod_form, x_min = NULL, n_indiv = 1, ll = 1.1, ul = 316, pars_to_get, scaling_var = 'dbh') {
  require(purrr)
  
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) 
  
  # Here, get the fitted log slope for all the predicted values. (x is dbh_pred and y is each column of fitted value)
  # Uses formula x/y dy/dx for log slope.
  log_slope <- function(x, y) (x[-1]/y[-1]) * diff(y)/diff(x)  
  
  if (variable == 'density' | variable == 'total_production') {
    pars_dens <- as.data.frame(do.call('cbind', extract(fit[['density']], pars_to_get[['density']])))
    
    if (dens_form == '1') {
      dens_fitted <- sapply(dbh_pred, pdf_pareto, xmin = x_min, alpha = pars_dens[,'alpha']) 
    }
    if (dens_form == '2') {
      dens_fitted <- sapply(dbh_pred, pdf_2part, xmin = x_min, alpha_low = pars_dens[,'alpha_low'], alpha_high = pars_dens[,'alpha_high'], tau = pars_dens[,'tau'])
    }
    if (dens_form == '3') {
      dens_fitted <- sapply(dbh_pred, pdf_3part, xmin = x_min, alpha_low = pars_dens[,'alpha_low'], alpha_mid = pars_dens[,'alpha_mid'], alpha_high = pars_dens[,'alpha_high'], tau_low = pars_dens[,'tau_low'], tau_high = pars_dens[,'tau_high'])
    }
    if (dens_form == 'w') {
      # Must manually rescale and remove upper and lower truncations
      trunc_pts <- pmap(pars_dens, function(m, n, ...) pweibull(q = c(ll,ul), shape = m, scale = n))
      dens_fitted <- sapply(dbh_pred, dweibull, shape = pars_dens[,'m'], scale = pars_dens[,'n']) 
      dens_fitted <- t(sapply(1:nrow(dens_fitted), function(i) {
        x <- dens_fitted[i,]/diff(trunc_pts[[i]])
        x
      }))
    }
    if (dens_form == 'ln') {
      dens_fitted <- sapply(dbh_pred, dlnorm, meanlog = pars_dens[,'mu_logn'], sdlog = pars_dens[,'sigma_logn']) 
    }
    dens_fitted <- dens_fitted * n_indiv
    dens_fittedslopes <- map_dfr(as.data.frame(t(dens_fitted)), ~ log_slope(dbh_pred, .))
    dens_fittedslopes_quant <- apply(dens_fittedslopes, 1, quantile, probs = qprobs, na.rm = TRUE)
  }
  
  if (variable == 'production' | variable == 'total_production') {
    pars_prod <- as.data.frame(do.call('cbind', extract(fit[['production']], pars_to_get[['production']])))
    if (prod_form == '1') {
      prod_fitted <- sapply(dbh_pred, powerlaw_log, beta0 = pars_prod[,'beta0'], beta1 = pars_prod[,'beta1'])
    }
    if (prod_form == '2') {
      prod_fitted <- sapply(dbh_pred, powerlaw_hinge_log, beta0 = pars_prod[,'beta0'], beta1_low = pars_prod[,'beta1_low'], beta1_high = pars_prod[,'beta1_high'], x0 = pars_prod[,'x0'], delta = pars_prod[,'delta'])
    }
    prod_fittedslopes <- map_dfr(as.data.frame(t(prod_fitted)), ~ log_slope(dbh_pred, .))
    prod_fittedslopes_quant <- apply(prod_fittedslopes, 1, quantile, probs = qprobs, na.rm = TRUE)
  }
  
  if (variable == 'total_production') {
    totalprod_fitted <- dens_fitted * prod_fitted
    
   
    totalprod_fittedslopes <- map_dfr(as.data.frame(t(totalprod_fitted)), ~ log_slope(dbh_pred, .))
    totalprod_fittedslopes_quant <- apply(totalprod_fittedslopes, 1, quantile, probs = qprobs, na.rm = TRUE)
  }
  
  if (variable == 'density') {
    out <- data.frame(dbh = dbh_pred[-1],
                      variable = rep('density', length(dbh_pred)-1),
                      t(dens_fittedslopes_quant))
  } 
  if (variable == 'production') {
    out <- data.frame(dbh = dbh_pred[-1],
                      variable = rep('production', length(dbh_pred)-1),
                      t(prod_fittedslopes_quant))
  }
  if (variable == 'total_production') {
    out <- data.frame(dbh = dbh_pred[-1],
                      variable = rep('total_production', length(dbh_pred)-1),
                      t(totalprod_fittedslopes_quant))
  }
  
  setNames(out, c(scaling_var, 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
}


# 4. Function to get the Bayesian R-squared and its quantiles (and bias correction factor)
# ----------------------------------------------------------------------------------------

bayesian_rsquared_production <- function(fit, x, y, prod_model) {
  # 3. Extract parameter estimates.
  production_par <- list('1' = c('beta0', 'beta1'),
						 '2' = c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta'))
  
  pars_to_get <- production_par[[prod_model]] 
  
  pars <- extract(fit, pars_to_get)
  pars <- as.data.frame(do.call(cbind, pars))
  
  # 4. Plug in dbh (x) to get posterior estimates of linear predictor of production
    
  # Take the log of the fitted values
  if (prod_model == '1') {
    prod_fitted <- log(do.call(rbind, pmap(pars, powerlaw_log, x = x)))
  } else {
    prod_fitted <- log(do.call(rbind, pmap(pars, powerlaw_hinge_log, x = x)))
  }

  # 5. Get residuals by subtracting log y from linear predictor
  resids <- -1 * sweep(prod_fitted, 2, log(y))
  
  # 6. Calculate variances and ratio
  pred_var <- apply(prod_fitted, 1, var)
  resid_var <- apply(resids, 1, var)
  r2s <- pred_var / (pred_var + resid_var)
  
  # 7. Bias correction factor
  # Sum of squared residuals
  ssq_resid <- apply(resids^2, 1, sum)
  # Standard error of estimates
  sse <- (ssq_resid / (length(y) - length(production_par[[prod_model]])))^0.5
  # Correction factors
  cfs <- exp((sse^2)/2)
  
  # Quantiles of rsq
  r2_quant <- quantile(r2s, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), na.rm = TRUE)
  r2_quant <- setNames(r2_quant, c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
  # Quantiles of correction factor
  cf_quant <- quantile(cfs, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), na.rm = TRUE)
  cf_quant <- setNames(cf_quant, c('q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975'))
  
  return(list(r2 = r2_quant, cf = cf_quant))
}

# 5. Master functions for the fits where density and production are done in separate models
# -----------------------------------------------------------------------------------------

# Just density: credible intervals of parameters, fitted values, fitted log slopes, WAIC, LOOIC
# Just production: credible intervals of parameters, fitted values, fitted log slopes, WAIC, LOOIC, Bayesian R-squared
# Density and production_combined: fitted values, fitted log slopes

extract_density <- function(dens_model, fg, year, xmin, n, use_subset = FALSE, n_chains = 3, scaling.var = 'dbh', fp = '~/forestlight/stanoutput', densityfitprefix = 'fit_density', scalingtype = 'production', LL = 1.1, UL = 316) {

  require(rstan)
  require(loo)
  require(Brobdingnag)

  # Load CSVs as stanfit object
  print('Loading stan fit . . .')
  files <- paste0(densityfitprefix, dens_model, '_', scalingtype, '_', fg, '_', year, '_', 1:n_chains, '.csv')
  fit <- list(density = read_stan_csv(file.path(fp,files)))
  
  # Get credible intervals of parameters
  print('Calculating credible intervals of parameters . . .')
  density_par <- list('1' = c('alpha'),
					  '2' = c('alpha_low', 'alpha_high', 'tau'),
					  '3' = c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high'),
					  'w' = c('m','n'),
					  'ln' = c('mu_logn', 'sigma_logn'))
  
  get_pars <- list(density = density_par[[dens_model]])
  
  param_cis <- cbind(data.frame(year = year, variable = 'density', model = dens_model, fg = fg), param_values(fit[['density']], get_pars[['density']]))
  
  # Get fitted values with credible intervals and prediction intervals
  print('Calculating fitted values and prediction intervals . . .')
  pred_interval <- fitted_predicted_values(fit, variable = 'density', dbh_pred, dens_form = dens_model, x_min = xmin, n_indiv = n, pars_to_get = get_pars, scaling_var = scaling.var, ll = LL, ul = UL)
  pred_interval <- data.frame(year = year, dens_model = dens_model, fg = fg, pred_interval)
  
  # Get fitted log slopes
  print('Calculating log slopes . . .')
  fitted_slopes <- fitted_slope_ci(fit, variable = 'density', dbh_pred, dens_form = dens_model, x_min = xmin, n_indiv = n, scaling_var = scaling.var, pars_to_get = get_pars, ll = LL, ul = UL)
  fitted_slopes <- data.frame(year = year, dens_model = dens_model, fg = fg, fitted_slopes)

  # Get WAIC and LOOIC
  ll_dens <- extract_log_lik(fit[['density']], 'log_lik')
  print('Calculating WAIC . . .')
  waic_dens <- waic(ll_dens)$estimates
  print('Calculating LOOIC . . .')
  loo_dens <- loo(ll_dens)$estimates

  list(waic = waic_dens, 
       loo = loo_dens, 
       param_cis = param_cis, 
       pred_interval = pred_interval,
	   fitted_slopes = fitted_slopes)
}

extract_production <- function(prod_model, fg, year, xmin, n, use_subset = FALSE, n_chains = 3, scaling.var = 'dbh', fp = '~/forestlight/stanoutput', fpdump = '~/forestlight/stanrdump', productionfitprefix = 'fit_production', scalingtype = 'production', dumpprefix = 'dump_', LL = 1.1, UL = 316) {
  require(rstan)
  require(loo)
  require(Brobdingnag)

  # Load CSVs as stanfit object
  print('Loading stan fit . . .')
  files <- paste0(productionfitprefix, prod_model, '_', scalingtype, '_', fg, '_', year, '_', 1:n_chains, '.csv')
  fit <- list(production = read_stan_csv(file.path(fp,files)))
  
  # Get credible intervals of parameters
  # Only get density parameters if production model is 1, only get production parameters if density model is 1
  print('Calculating credible intervals of parameters . . .')
  production_par <- list('1' = c('beta0', 'beta1', 'sigma'),
						 '2' = c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta', 'sigma'))
  
  get_pars <- list(production = production_par[[prod_model]])
  
  param_cis <- cbind(data.frame(year = year, variable = 'production', model = prod_model, fg = fg), param_values(fit[['production']], get_pars[['production']]))

  # Get fitted values with credible intervals and prediction intervals
  print('Calculating fitted values and prediction intervals . . .')
  pred_interval <- fitted_predicted_values(fit, variable = 'production', dbh_pred, prod_form = prod_model, x_min = xmin, n_indiv = n, pars_to_get = get_pars, scaling_var = scaling.var, ll = LL, ul = UL)
  pred_interval <- data.frame(year = year, prod_model = prod_model, fg = fg, pred_interval)
  
  # Get fitted log slopes
  print('Calculating log slopes . . .')
  fitted_slopes <- fitted_slope_ci(fit, variable = 'production', dbh_pred, prod_form = prod_model, x_min = xmin, n_indiv = n, scaling_var = scaling.var, pars_to_get = get_pars, ll = LL, ul = UL)
  fitted_slopes <- data.frame(year = year, prod_model = prod_model, fg = fg, fitted_slopes)

  # Get WAIC and LOOIC
  ll_prod <- extract_log_lik(fit[['production']], 'log_lik')
  print('Calculating WAIC . . .')
  waic_prod <- waic(ll_prod)$estimates
  print('Calculating LOOIC . . .')
  loo_prod <- loo(ll_prod)$estimates
  
  # Calculate R-squared and correction factor
  print('Calculating Bayesian R-squared and bias CF . . .')
  # Load the dump file from the model so that the R2 can be calculated
  dumpfile <- paste0(dumpprefix, fg, '_', year, '.r')

  source(file.path(fpdump, dumpfile)) # Creates variables x and y.
  
  r2s_cfs <- bayesian_rsquared_production(fit[['production']], x, y, prod_model)
  
  list(waic = waic_prod, 
       loo = loo_prod, 
       param_cis = param_cis, 
       pred_interval = pred_interval,
	   fitted_slopes = fitted_slopes,
	   r2s = r2s_cfs$r2,
	   cfs = r2s_cfs$cf)
}

extract_totalproduction <- function(dens_model, prod_model, fg, year, xmin, n, use_subset = FALSE, n_chains = 3, scaling.var = 'dbh', fp = '~/forestlight/stanoutput', densityfitprefix = 'fit_density', productionfitprefix = 'fit_production', scalingtype = 'production', LL = 1.1, UL = 316) {
  require(rstan)
  require(Brobdingnag)
  
  # Load CSVs as stanfit object
  print('Loading stan fit . . .')
  files_density <- paste0(densityfitprefix, dens_model, '_', 'production', '_', fg, '_', year, '_', 1:n_chains, '.csv')
  files_production <- paste0(productionfitprefix, prod_model, '_', scalingtype, '_', fg, '_', year, '_', 1:n_chains, '.csv')
  fit <- list(density = read_stan_csv(file.path(fp,files_density)), production = read_stan_csv(file.path(fp,files_production)))
  
  # Get fitted values with credible intervals and prediction intervals
  print('Calculating fitted values and prediction intervals . . .')
  
  density_par <- list('1' = c('alpha'),
                      '2' = c('alpha_low', 'alpha_high', 'tau'),
                      '3' = c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high'),
                      'w' = c('m','n'),
                      'ln' = c('mu_logn', 'sigma_logn'))
  production_par <- list('1' = c('beta0', 'beta1', 'sigma'),
                         '2' = c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta', 'sigma'))					  
  
  get_pars <- list(density = density_par[[dens_model]], production = production_par[[prod_model]])
  
  pred_interval <- fitted_predicted_values(fit, variable = 'total_production', dbh_pred, dens_form = dens_model, prod_form = prod_model, x_min = xmin, n_indiv = n, pars_to_get = get_pars, scaling_var = scaling.var, ll = LL, ul = UL)
  pred_interval <- data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, pred_interval)
  
  # Get fitted log slopes
  print('Calculating log slopes . . .')
  fitted_slopes <- fitted_slope_ci(fit, variable = 'total_production', dbh_pred, dens_form = dens_model, prod_form = prod_model, x_min = xmin, n_indiv = n, scaling_var = scaling.var, pars_to_get = get_pars, ll = LL, ul = UL)
  fitted_slopes <- data.frame(year = year, dens_model = dens_model, prod_model = prod_model, fg = fg, fitted_slopes)
  
  list(pred_interval = pred_interval,
       fitted_slopes = fitted_slopes)
}

