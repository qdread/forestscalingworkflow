# Model output extraction functions for richness

bayesian_rsquared_richness <- function(fit, x, y, model) {
  # 3. Extract parameter estimates.
   <- list('linear' = c('beta0', 'beta1'),
						 'mixed' = c('beta0', 'beta1_low', 'beta1_high', 'x0', 'delta'))
  
  if ()
  
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