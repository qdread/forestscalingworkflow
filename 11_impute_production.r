# Imputed values for new recruit individuals that did not have production measured because they only had one diameter measurement
# Must be run BEFORE creating log bins and other visualization data
# QDR / Forestlight / 13 Jan 2020


# Load data ---------------------------------------------------------------

library(tidyverse)
library(forestscaling)

# Load raw data and the parameter values 

load(file.path(gdrive_path, 'data/rawdataobj_alternativecluster.r'))

params <- read_csv(file.path(gdrive_path, 'finalcsvs/piecewise_paramci_by_fg.csv')) # parameter values
cfs <- read_csv(file.path(gdrive_path, 'finalcsvs/piecewise_cf_by_fg.csv')) # bias correction factors

# Reshape parameter and correction factor DFs into wide format.
prodpars_wide <- params %>%
  filter(variable == 'production', !fg %in% 'unclassified', !parameter %in% 'sigma') %>%
  mutate(parameter = if_else(parameter == 'beta0' & model == 2, 'beta0_2segment', parameter)) %>%
  select(fg, parameter, q50) %>%
  pivot_wider(names_from = parameter, values_from = q50)

cfs_wide <- cfs %>%
  filter(!fg %in% 'unclassified') %>%
  select(prod_model, fg, q50) %>%
  pivot_wider(names_from = prod_model, values_from = q50, names_prefix = 'cf') 

pars_cfs <- left_join(prodpars_wide, cfs_wide)


# Define function ---------------------------------------------------------

# For all rows that are new recruit, calculate the imputed production value for 1 and 2 segment production models.

# We will use group-specific regression parameters, median values.
# For unclassified trees, use the all-tree regression.

# define function to be used on each year's data
# version that only returns the two columns
impute_production <- function(dat) {
  # Join data with parameters and correction factors
  dat <- dat %>% 
    select(fg, dbh_corr, recruit) %>%
    mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'alltree')) %>%
    left_join(pars_cfs)
  # Apply the two functions and multiply by correction factor to all the dbhs (including non-recruits)
  dat %>%
    transmute(prod_imputed_1segment = powerlaw_log(x = dbh_corr, beta0 = beta0, beta1 = beta1) * cf1,
              prod_imputed_2segment = powerlaw_hinge_log(x = dbh_corr, x0 = x0, beta0 = beta0_2segment, beta1_low = beta1_low, beta1_high = beta1_high, delta = delta) * cf2)
  
}

# version that returns a single column, to be used with mutate() 
impute_production_bymodel <- function(dat, mod) {
  # Join data with parameters and correction factors
  dat <- dat %>% 
    select(fg, dbh_corr, production, recruit) %>%
    mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'alltree')) %>%
    left_join(pars_cfs)
  # Apply the two functions and multiply by correction factor to all the dbhs (including non-recruits)
  dat <- dat %>%
    mutate(prod_imputed_1segment = powerlaw_log(x = dbh_corr, beta0 = beta0, beta1 = beta1) * cf1,
              prod_imputed_2segment = powerlaw_hinge_log(x = dbh_corr, x0 = x0, beta0 = beta0_2segment, beta1_low = beta1_low, beta1_high = beta1_high, delta = delta) * cf2) %>%
    mutate(prod_imputed_1segment = if_else(recruit, prod_imputed_1segment, production),
           prod_imputed_2segment = if_else(recruit, prod_imputed_2segment, production))
  if (mod == 1) {
    return(dat$prod_imputed_1segment)
  } else {
    return(dat$prod_imputed_2segment)
  }
}

# Apply function to all DFs -----------------------------------------------

alltree_light_90 <- alltree_light_90 %>%
  mutate(production_imputed1 = impute_production_bymodel(., 1),
         production_imputed2 = impute_production_bymodel(., 2))
alltree_light_95 <- alltree_light_95 %>%
  mutate(production_imputed1 = impute_production_bymodel(., 1),
         production_imputed2 = impute_production_bymodel(., 2))

alltreedat <- map(alltreedat, function(x) x %>% mutate(production_imputed1 = impute_production_bymodel(x, 1),
                                                       production_imputed2 = impute_production_bymodel(x, 2)))

light_fg_90 <- map(light_fg_90, function(x) x %>% mutate(production_imputed1 = impute_production_bymodel(x, 1),
                                                         production_imputed2 = impute_production_bymodel(x, 2)))

light_fg_95 <- map(light_fg_95, function(x) x %>% mutate(production_imputed1 = impute_production_bymodel(x, 1),
                                                         production_imputed2 = impute_production_bymodel(x, 2)))

fgdat <- map(fgdat, function(y) map(y, function(x) x %>% mutate(production_imputed1 = impute_production_bymodel(x, 1),
                                                                production_imputed2 = impute_production_bymodel(x, 2))))

save(alltree_light_90, alltree_light_95, alltreedat, light_fg_90, light_fg_95, fgdat, file = file.path(gdrive_path, 'data/rawdataobj_withimputedproduction.RData'))
