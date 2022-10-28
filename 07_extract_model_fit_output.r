# For each CMDstan model fit, get the following
# Credible intervals of parameters
# Credible and prediction intervals to make graphs
# Information criteria
# Fitted slopes in log space
# Bayesian R-squared

# Edited 04 October 2022 to remove fits not present in current version of manuscript.

# DENSITY - PRODUCTION SCALINGS
# =============================

source('R_functions/model_output_extraction_functions.r')

library(purrr)
library(dplyr)
library(foreach)
library(doParallel)

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
					   					   
mod_df <- expand.grid(variable = 'total_production',
					  dens_model = 1:3,
                      prod_model = 1:2,
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = 1995, 
                      stringsAsFactors = FALSE)

mod_df <- rbind(dens_df, prod_df, mod_df)
					  
min_n <- read.csv('stanrdump/min_n.csv', stringsAsFactors = FALSE)

mod_df <- mod_df %>%
  left_join(min_n)

dbh_pred <- exp(seq(log(1), log(315), length.out = 101))


registerDoParallel(cores = 8)

tmp <- foreach(i = 1:nrow(mod_df)) %dopar% {

	if (mod_df$variable[i] == 'density') {
		fit_info <- extract_density(dens_model = mod_df$dens_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									use_subset = FALSE)
	}
	if (mod_df$variable[i] == 'production') {
		fit_info <- extract_production(prod_model = mod_df$prod_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									dumpprefix = 'dump_production_',
									use_subset = FALSE)
	}
	if (mod_df$variable[i] == 'total_production') {
		fit_info <- extract_totalproduction(dens_model = mod_df$dens_model[i],
									prod_model = mod_df$prod_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									use_subset = FALSE)
	}

	save(fit_info, file = paste0('stanoutput/fitinfo/pw_info_',mod_df$variable[i],'_',i,'.r'))
	message('Fit ', i, ' saved')

}

# GROWTH AS DIAMETER PER TIME SCALING
# INDIVIDUAL PRODUCTION ONLY
# ===================================

source('R_functions/model_output_extraction_functions.r')

library(purrr)
library(dplyr)
library(foreach)
library(doParallel)

prod_df <- expand.grid(variable = 'production',
					   dens_model = as.numeric(NA),
					   prod_model = 1:2,
					   fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
					   year = 1995,
					   stringsAsFactors = FALSE) 
					   
min_n <- read.csv('stanrdump/min_n.csv', stringsAsFactors = FALSE)	


mod_df <- prod_df %>%
  left_join(min_n)

dbh_pred <- exp(seq(log(1), log(315), length.out = 101))

registerDoParallel(cores = 8)

tmp <- foreach(i = 1:nrow(mod_df)) %dopar% {
		fit_info <- extract_production(prod_model = mod_df$prod_model[i],
									fg = mod_df$fg[i],
									year = mod_df$year[i],
									xmin = mod_df$xmin[i],
									n = mod_df$n[i],
									scalingtype = 'diamgrowthscaling',
									dumpprefix = 'dump_diamgrowthscaling_',
									use_subset = FALSE)
	
	save(fit_info, file = paste0('stanoutput/fitinfo/diamgrowthpw_info_',i,'.r'))
	message('Fit ', i, ' saved')

}				   
