# Load all the fit info into a big list and put each thing into its own data frame
# (fit info from the cmdstan fits)
# Version created for piecewise. Only done 1995.

# Alternate version created for the ones where density and production are fit separately.
# Edited 17 June 2019 to accommodate the ones with total light and crown volume used as scaling.

# DENSITY-PRODUCTION SCALINGS
# ===========================

fp <- 'stanoutput/fitinfo'

library(dplyr)
library(purrr)

dens_df <- expand.grid(variable = 'density',
					   dens_model = 1:3,
					   prod_model = as.numeric(NA),
					   fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
					   year = 1995,
					   stringsAsFactors = FALSE)

prod_df <- expand.grid(variable = 'production',
					   dens_model = as.numeric(NA),
					   prod_model = 1:2,
					   fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
					   year = 1995,
					   stringsAsFactors = FALSE)
					   					   
totalprod_df <- expand.grid(variable = 'total_production',
					  dens_model = 1:3,
                      prod_model = 1:2,
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = 1995, 
                      stringsAsFactors = FALSE)

mod_df <- rbind(dens_df, prod_df, totalprod_df)
					  

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(file.path(fp, paste0('pw_info_',mod_df$variable[i],'_',i,'.r')))
  fit_info
})

# Extract WAIC and LOOIC for density and production from each one.
# ----------------------------------------------------------------

get_ics <- function(x) {
	ics <- data.frame(criterion = c('WAIC', 'LOOIC'),
					  IC_value = c(x$waic['waic','Estimate'], x$loo['looic','Estimate']),
					  IC_stderr = c(x$waic['waic','SE'], x$loo['looic','SE']))
	return(ics)
}

idx <- which(mod_df$variable %in% c('density', 'production'))
fit_ics <- map_dfr(fit_info_list[idx], get_ics)
fit_ics <- cbind(mod_df[idx,][rep(1:nrow(mod_df[idx,]), each=2),], fit_ics)
write.csv(fit_ics, file = 'finalcsvs/piecewise_ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- map_dfr(fit_info_list[idx], 'param_cis')

write.csv(param_cis, file = 'finalcsvs/piecewise_paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- map_dfr(fit_info_list, 'pred_interval')

write.csv(pred_values, file = 'finalcsvs/piecewise_ci_by_fg.csv', row.names = FALSE)

# Combine fitted slopes into single data frame.
# ---------------------------------------------

fitted_slopes <- map_dfr(fit_info_list, 'fitted_slopes')

write.csv(fitted_slopes, file = 'finalcsvs/piecewise_fitted_slopes_by_fg.csv', row.names = FALSE)

# Combine R-squared values into single data frame.
# ------------------------------------------------

idx_p <- which(mod_df$variable %in% c('production'))
r2s <- do.call(rbind, map(fit_info_list[idx_p], 'r2s'))
r2s <- cbind(prod_df, r2s)

write.csv(r2s, file = 'finalcsvs/piecewise_r2_by_fg.csv', row.names = FALSE)

# Combine bias correction factors into single data frame.
# -------------------------------------------------------
 
cfs <- do.call(rbind, map(fit_info_list[idx_p], 'cfs'))
cfs <- cbind(prod_df, cfs)

write.csv(cfs, file = 'finalcsvs/piecewise_cf_by_fg.csv', row.names = FALSE)


# DENSITY-TOTAL LIGHT SCALINGS
# ============================

mod_df <- rbind(prod_df, totalprod_df)
					  

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(file.path(fp, paste0('lightpw_info_',mod_df$variable[i],'_',i,'.r')))
  fit_info
})

# Extract WAIC and LOOIC for total light from each one.
# ----------------------------------------------------------------

get_ics <- function(x) {
	ics <- data.frame(criterion = c('WAIC', 'LOOIC'),
					  IC_value = c(x$waic['waic','Estimate'], x$loo['looic','Estimate']),
					  IC_stderr = c(x$waic['waic','SE'], x$loo['looic','SE']))
	return(ics)
}

idx <- which(mod_df$variable %in% c('production'))
fit_ics <- map_dfr(fit_info_list[idx], get_ics)
fit_ics <- cbind(mod_df[idx,][rep(1:nrow(mod_df[idx,]), each=2),], fit_ics)
fit_ics$variable <- 'incoming light individual'
write.csv(fit_ics, file = 'finalcsvs/light_piecewise_ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- map_dfr(fit_info_list[idx], 'param_cis')
param_cis$variable <- 'incoming light individual'

write.csv(param_cis, file = 'finalcsvs/light_piecewise_paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- map_dfr(fit_info_list, 'pred_interval')
pred_values$variable <- factor(pred_values$variable, levels = c('production', 'production_fitted', 'total_production', 'total_production_fitted'), labels = c('incoming_light', 'incoming_light_fitted', 'total_incoming_light', 'total_incoming_light_fitted'))

write.csv(pred_values, file = 'finalcsvs/light_piecewise_ci_by_fg.csv', row.names = FALSE)

# Combine fitted slopes into single data frame.
# ---------------------------------------------

fitted_slopes <- map_dfr(fit_info_list, 'fitted_slopes')
fitted_slopes$variable <- factor(fitted_slopes$variable, levels = c('production', 'total_production'), labels = c('incoming_light', 'total_incoming_light'))

write.csv(fitted_slopes, file = 'finalcsvs/light_piecewise_fitted_slopes_by_fg.csv', row.names = FALSE)

# Combine R-squared values into single data frame.
# ------------------------------------------------

idx_p <- which(mod_df$variable %in% c('production'))
r2s <- do.call(rbind, map(fit_info_list[idx_p], 'r2s'))
r2s <- cbind(prod_df, r2s)
r2s$variable <- 'incoming light individual'

write.csv(r2s, file = 'finalcsvs/light_piecewise_r2_by_fg.csv', row.names = FALSE)

# Combine bias correction factors into single data frame.
# -------------------------------------------------------
 
cfs <- do.call(rbind, map(fit_info_list[idx_p], 'cfs'))
cfs <- cbind(prod_df, cfs)
cfs$variable <- 'incoming light individual'

write.csv(cfs, file = 'finalcsvs/light_piecewise_cf_by_fg.csv', row.names = FALSE)

# DENSITY x CROWN VOLUME SCALING
# ==============================

mod_df <- rbind(prod_df, totalprod_df)
					  

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(file.path(fp, paste0('volumepw_info_',mod_df$variable[i],'_',i,'.r')))
  fit_info
})

# Extract WAIC and LOOIC for total light from each one.
# ----------------------------------------------------------------

get_ics <- function(x) {
	ics <- data.frame(criterion = c('WAIC', 'LOOIC'),
					  IC_value = c(x$waic['waic','Estimate'], x$loo['looic','Estimate']),
					  IC_stderr = c(x$waic['waic','SE'], x$loo['looic','SE']))
	return(ics)
}

idx <- which(mod_df$variable %in% c('production'))
fit_ics <- map_dfr(fit_info_list[idx], get_ics)
fit_ics <- cbind(mod_df[idx,][rep(1:nrow(mod_df[idx,]), each=2),], fit_ics)
fit_ics$variable <- 'crown volume individual'
write.csv(fit_ics, file = 'finalcsvs/volume_piecewise_ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- map_dfr(fit_info_list[idx], 'param_cis')
param_cis$variable <- 'crown volume individual'

write.csv(param_cis, file = 'finalcsvs/volume_piecewise_paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- map_dfr(fit_info_list, 'pred_interval')
pred_values$variable <- factor(pred_values$variable, levels = c('production', 'production_fitted', 'total_production', 'total_production_fitted'), labels = c('crown_volume', 'crown_volume_fitted', 'total_crown_volume', 'total_crown_volume_fitted'))

write.csv(pred_values, file = 'finalcsvs/volume_piecewise_ci_by_fg.csv', row.names = FALSE)

# Combine fitted slopes into single data frame.
# ---------------------------------------------

fitted_slopes <- map_dfr(fit_info_list, 'fitted_slopes')
fitted_slopes$variable <- factor(fitted_slopes$variable, levels = c('production', 'total_production'), labels = c('crown_volume', 'total_crown_volume'))

write.csv(fitted_slopes, file = 'finalcsvs/volume_piecewise_fitted_slopes_by_fg.csv', row.names = FALSE)

# Combine R-squared values into single data frame.
# ------------------------------------------------

idx_p <- which(mod_df$variable %in% c('production'))
r2s <- do.call(rbind, map(fit_info_list[idx_p], 'r2s'))
r2s <- cbind(prod_df, r2s)
r2s$variable <- 'crown volume individual'

write.csv(r2s, file = 'finalcsvs/volume_piecewise_r2_by_fg.csv', row.names = FALSE)

# Combine bias correction factors into single data frame.
# -------------------------------------------------------
 
cfs <- do.call(rbind, map(fit_info_list[idx_p], 'cfs'))
cfs <- cbind(prod_df, cfs)
cfs$variable <- 'crown volume individual'

write.csv(cfs, file = 'finalcsvs/volume_piecewise_cf_by_fg.csv', row.names = FALSE)


# GROWTH AS DIAMETER PER TIME SCALINGS
# ADDED 24 JULY 2019
# ====================================

mod_df <- prod_df

fit_info_list <- map(1:nrow(mod_df), function(i) {
  load(file.path(fp, paste0('diamgrowthpw_info_',i,'.r')))
  fit_info
})

# Extract WAIC and LOOIC for diameter growth from each one.
# ----------------------------------------------------------------

get_ics <- function(x) {
	ics <- data.frame(criterion = c('WAIC', 'LOOIC'),
					  IC_value = c(x$waic['waic','Estimate'], x$loo['looic','Estimate']),
					  IC_stderr = c(x$waic['waic','SE'], x$loo['looic','SE']))
	return(ics)
}

fit_ics <- map_dfr(fit_info_list, get_ics)
fit_ics <- cbind(mod_df[rep(1:nrow(mod_df), each=2),], fit_ics)
fit_ics$variable <- 'diameter growth individual'
write.csv(fit_ics, file = 'finalcsvs/diamgrowth_piecewise_ics_by_fg.csv', row.names = FALSE)

# Combine parameter credible intervals into single data frame.
# ------------------------------------------------------------

param_cis <- map_dfr(fit_info_list, 'param_cis')
param_cis$variable <- 'diameter growth individual'

write.csv(param_cis, file = 'finalcsvs/diamgrowth_piecewise_paramci_by_fg.csv', row.names = FALSE)

# Combine predicted values into single data frame.
# ------------------------------------------------

pred_values <- map_dfr(fit_info_list, 'pred_interval')
pred_values$variable <- factor(pred_values$variable, levels = c('production', 'production_fitted'), labels = c('diameter_growth', 'diameter_growth_fitted'))

write.csv(pred_values, file = 'finalcsvs/diamgrowth_piecewise_ci_by_fg.csv', row.names = FALSE)

# Combine fitted slopes into single data frame.
# ---------------------------------------------

fitted_slopes <- map_dfr(fit_info_list, 'fitted_slopes')
fitted_slopes$variable <- 'diameter growth individual'

write.csv(fitted_slopes, file = 'finalcsvs/diamgrowth_piecewise_fitted_slopes_by_fg.csv', row.names = FALSE)

# Combine R-squared values into single data frame.
# ------------------------------------------------

r2s <- do.call(rbind, map(fit_info_list, 'r2s'))
r2s <- cbind(prod_df, r2s)
r2s$variable <- 'diameter growth individual'

write.csv(r2s, file = 'finalcsvs/diamgrowth_piecewise_r2_by_fg.csv', row.names = FALSE)

# Combine bias correction factors into single data frame.
# -------------------------------------------------------
 
cfs <- do.call(rbind, map(fit_info_list, 'cfs'))
cfs <- cbind(prod_df, cfs)
cfs$variable <- 'diameter growth individual'

write.csv(cfs, file = 'finalcsvs/diamgrowth_piecewise_cf_by_fg.csv', row.names = FALSE)

