# Get credible intervals for LUE slopes, and fitted values for the LUE trend (for plotting).
# Use only 1 segment models.
# QDR / Forestlight / 03 May 2021

get_lue_info <- function(fg, year = 1995, n_chains = 3) {
  require(rstan)
  require(Brobdingnag)
  require(purrr)
  require(forestscaling)
  
  qprobs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) 

  # names of parameters
  production_par <- c('beta0', 'beta1', 'sigma')			  
  
  message(paste("Functional group", fg))
  message('Loading stan fit . . .')
  
  files_lightcaptured <- glue::glue('~/forestlight/stanoutput/fit_production1_lightcapturedscaling_{fg}_{year}_{1:n_chains}.csv')
  files_production <- glue::glue('~/forestlight/stanoutput/fit_production1_production_{fg}_{year}_{1:n_chains}.csv')
  fit <- list(lightcaptured = read_stan_csv(files_lightcaptured), 
              production = read_stan_csv(files_production))
  
  message('Extracting parameters and calculating CIs . . .')

  pars_lightcaptured <- as.data.frame(do.call('cbind', extract(fit[['lightcaptured']], production_par)))
  pars_production <- as.data.frame(do.call('cbind', extract(fit[['production']], production_par)))
  
  # Get difference of slope and its quantiles
  lue_slope <- pars_production$beta1 - pars_lightcaptured$beta1
  lue_slope_quant <- quantile(lue_slope, probs = qprobs, na.rm = TRUE)

  message('Getting fitted values . . .')

  lightcaptured_fitted <- sapply(dbh_pred, powerlaw_log, beta0 = pars_lightcaptured[,'beta0'], beta1 = pars_lightcaptured[,'beta1'])
  prod_fitted <- sapply(dbh_pred, powerlaw_log, beta0 = pars_production[,'beta0'], beta1 = pars_production[,'beta1'])
  
  message('Loading raw data and calculating bias correction factor . . .')
  
  source(glue::glue('~/forestlight/stanrdump/dump_production_{fg}_{year}.r'))
  prod_raw <- list(x = x, y = y)
  source(glue::glue('~/forestlight/stanrdump/dump_lightcapturedscaling_{fg}_{year}.r'))
  lightcaptured_raw <- list(x = x, y = y)
  
  # We need all fitted values for production, not just the 101 values, to calculate correction factor.
  lightcaptured_fitted_all <- pmap(pars_lightcaptured[, c('beta0','beta1')], powerlaw_log, x = lightcaptured_raw$x)
  prod_fitted_all <- pmap(pars_production[, c('beta0','beta1')], powerlaw_log, x = prod_raw$x)
  
  lightcaptured_cf <- corr_factor(y = lightcaptured_raw$y, y_fit = lightcaptured_fitted_all, n_pars = 2)
  prod_cf <- corr_factor(y = prod_raw$y, y_fit = prod_fitted_all, n_pars = 2)
  
  lightcaptured_fitted <- sweep(lightcaptured_fitted, 1, lightcaptured_cf, `*`)
  prod_fitted <- sweep(prod_fitted, 1, prod_cf, `*`)
  
  message('Getting quantiles of fitted values (almost done!) . . .')

  # Quotient of fitted values is light use efficiency
  lue_fitted <- prod_fitted / lightcaptured_fitted
  
  lue_fitted_quant <- apply(lue_fitted, 2, quantile, probs = qprobs, na.rm = TRUE)
  
  lue_fitted_quant <- data.frame(dbh = dbh_pred,
                                 variable = 'light use efficiency',
                                 lue_fitted_quant)
  
  names(lue_fitted_quant) <- c('dbh', 'variable', 'q025', 'q05', 'q25', 'q50', 'q75', 'q95', 'q975')
  
  return(list(lue_slope_quant = lue_slope_quant, lue_fitted_quant = lue_fitted_quant))
  
}

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


# Execute function for each FG ----------------------

library(purrr)

fgs <- c('alltree', paste0('fg', 1:5))
dbh_pred <- exp(seq(log(1), log(315), length.out = 101))

lue_info <- map(fgs, get_lue_info)
save(lue_info, file = '~/forestlight/stanoutput/fitinfo/lue_info.RData')