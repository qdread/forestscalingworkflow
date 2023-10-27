# Create piecewise data for plotting.
# Modified 13 Jan 2020: Use imputed production totals (1 segment production) to create total production bins (but not individual)
# Modified 09 Jan 2020: Fix incorrect number of individuals in some of the observed data.
# Modified 13 Dec 2019: Create both observed and predicted data in this script. 
# Modified 13 Dec 2019: Use the new up to date correction factor based on the true correction of Jensen's Inequality.
# Modified 02 Jul 2019: Calculate new correction factors to get the proper normalizations (based on integral with bounded limits)

library(tidyverse)
library(forestscaling)

fp_out <- 'ForestLight/data/data_forplotting'

years <- seq(1990, 2010, 5)
fg_names <- c("fg1", "fg2", "fg3", "fg4", "fg5", "unclassified")

# Load all the by-year binned data.
load('ForestLight/data/data_binned/bin_object_singleyear.RData')

# Load raw data
load('ForestLight/data/rawdataobj_withimputedproduction.RData')

binedgedata <- densitybin_byyear %>% filter(fg == 'all', year == 1995) 
area_core <- 42.84

# Light received, crown volume, light captured, and leaf area bins (the latter two added 27 Apr 2021)
totallightbins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = light_received, edges = binedgedata))
totalvolbins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = crownvolume, edges = binedgedata))
totallightcapturedbins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = light_captured, edges = binedgedata))
totalleafareabins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = leaf_area, edges = binedgedata))

totallightbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$light_received, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totallightbins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)
totalvolbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$crownvolume, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totalvolbins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)
totallightcapturedbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$light_captured, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totallightcapturedbins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)
totalleafareabins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$leaf_area, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totalleafareabins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)

indivlightbins_all <- alltree_light_95 %>%
  mutate(indivlight_bin = cut(dbh_corr, breaks = c(binedgedata$bin_min[1], binedgedata$bin_max), labels = binedgedata$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(indivlight_bin) %>%
  do(c(n = nrow(.), quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
indivlightbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  mutate(indivlight_bin =cut(dbh_corr, breaks = c(binedgedata$bin_min[1], binedgedata$bin_max), labels = binedgedata$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, indivlight_bin) %>%
  do(c(n = nrow(.), quantile(.$light_received, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
indivlightbins_fg <- data.frame(fg = 'all', indivlightbins_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(indivlightbins_fg)) %>%
  mutate(indivlight_bin = as.numeric(as.character(indivlight_bin))) %>%
  rename(bin_midpoint = indivlight_bin)

indivlightcapturedbins_all <- alltree_light_95 %>%
  mutate(indivlightcaptured_bin = cut(dbh_corr, breaks = c(binedgedata$bin_min[1], binedgedata$bin_max), labels = binedgedata$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(indivlightcaptured_bin) %>%
  do(c(n = nrow(.), quantile(.$light_captured, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
indivlightcapturedbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  mutate(indivlightcaptured_bin =cut(dbh_corr, breaks = c(binedgedata$bin_min[1], binedgedata$bin_max), labels = binedgedata$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, indivlightcaptured_bin) %>%
  do(c(n = nrow(.), quantile(.$light_captured, c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('bin_count', 'q025','q25','q50','q75','q975')))
indivlightcapturedbins_fg <- data.frame(fg = 'all', indivlightcapturedbins_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(indivlightcapturedbins_fg)) %>%
  mutate(indivlightcaptured_bin = as.numeric(as.character(indivlightcaptured_bin))) %>%
  rename(bin_midpoint = indivlightcaptured_bin)

write.csv(indivlightbins_fg, file.path(fp_out, 'obs_indivlight.csv'), row.names = FALSE)
write.csv(totallightbins_fg, file.path(fp_out, 'obs_totallight.csv'), row.names = FALSE)
write.csv(indivlightcapturedbins_fg, file.path(fp_out, 'obs_indivlightcaptured.csv'), row.names = FALSE)
write.csv(totallightcapturedbins_fg, file.path(fp_out, 'obs_totallightcaptured.csv'), row.names = FALSE)
write.csv(totalvolbins_fg, file.path(fp_out, 'obs_totalvol.csv'), row.names = FALSE)
write.csv(totalleafareabins_fg, file.path(fp_out, 'obs_totalleafarea.csv'), row.names = FALSE)																									 
# Observed individual production
obs_indivprod_df <- map(alltreedat[-1],
                         function(x) with(x %>% filter(!recruit), cloudbin_across_years(dat_values = production, dat_classes = dbh_corr, edges = binedgedata, n_census = 1)))
obs_indivprod_df <- cbind(fg = 'all', year = rep(years, each = nrow(binedgedata)), do.call(rbind, obs_indivprod_df))

obs_indivprod_fg <- replicate(length(fg_names), list())

for (i in 1:length(fg_names)) {
  for (j in 1:length(years)) {
    obs_indivprod_fg[[i]][[j]] <-
      cbind(fg = fg_names[i], year = years[j],
            cloudbin_across_years(dat_values = fgdat[[i]][[j+1]] %>% filter(!recruit) %>% pull(production),
                                  dat_classes = fgdat[[i]][[j+1]] %>% filter(!recruit) %>% pull(dbh_corr),
                                  edges = binedgedata,
                                  n_census = 1))
  }
}

obs_indivprod_fg <- do.call(rbind, map(obs_indivprod_fg, function(x) do.call(rbind, x)))
obs_indivprod_df <- rbind(obs_indivprod_df, obs_indivprod_fg)

# Load correction factor for total production
cf_production <- read_csv('finalcsvs/piecewise_cf_by_fg.csv') %>% 
  select(prod_model, fg, q50) %>% 
  rename(corr_factor = q50) %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))


# Convert to individuals and production per hectare.
obs_dens <- densitybin_byyear %>%
  mutate(bin_value = bin_value / area_core)
obs_totalprod <- totalproductionbin_byyear %>%
  mutate(bin_value = bin_value / area_core)


ci_df <- read.csv('finalcsvs/piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

pred_dens <- ci_df %>%
  filter(variable == 'density') %>%
  select(-variable) %>%
  mutate_at(vars(starts_with('q')), ~ ./area_core)

fitted_indivprod <- ci_df %>%
  filter(variable == 'production_fitted') %>%
  select(-variable)

fitted_totalprod <- ci_df %>%
  filter(variable == 'total_production_fitted') %>%
  select(-variable) 

pred_indivprod <- ci_df %>%
  filter(variable == 'production') %>%
  select(-variable)

pred_totalprod <- ci_df %>%
  filter(variable == 'total_production') %>%
  select(-variable) 

fitted_totalprod <- fitted_totalprod %>%
  left_join(cf_production) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

pred_totalprod <- pred_totalprod %>%
  left_join(cf_production) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

# Added 09 Jan 2020: replace incorrect bin counts in observed production data with the bin counts from observed density.
joined_counts <- obs_indivprod_df %>% 
  select(fg, year, bin_midpoint) %>%
  left_join(obs_dens %>% select(fg, year, bin_midpoint, bin_count))

obs_indivprod_df$mean_n_individuals <- joined_counts$bin_count

write.csv(obs_dens, file.path(fp_out, 'obs_dens.csv'), row.names = FALSE)
write.csv(obs_totalprod, file.path(fp_out, 'obs_totalprod.csv'), row.names = FALSE)
write.csv(obs_indivprod_df, file.path(fp_out, 'obs_indivprod.csv'), row.names = FALSE)
write.csv(pred_dens, file.path(fp_out, 'pred_dens.csv'), row.names = FALSE)
write.csv(pred_indivprod, file.path(fp_out, 'pred_indivprod.csv'), row.names = FALSE)
write.csv(pred_totalprod, file.path(fp_out, 'pred_totalprod.csv'), row.names = FALSE)
write.csv(fitted_indivprod, file.path(fp_out, 'fitted_indivprod.csv'), row.names = FALSE)
write.csv(fitted_totalprod, file.path(fp_out, 'fitted_totalprod.csv'), row.names = FALSE)

#################
# For individual + total light scaling and total volume scaling
# Added 20 June 2019

library(dplyr)

ci_df <- read.csv('finalcsvs/light_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

cf_totallight<- read_csv('finalcsvs/light_piecewise_cf_by_fg.csv') %>% 
  select(prod_model, fg, q50) %>% 
  rename(corr_factor = q50) %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

fitted_indivlight <- ci_df %>%
  filter(variable == 'incoming_light_fitted') %>%
  select(-variable)

fitted_totallight <- ci_df %>%
  filter(variable == 'total_incoming_light_fitted') %>%
  select(-variable) 

pred_indivlight <- ci_df %>%
  filter(variable == 'incoming_light') %>%
  select(-variable)

pred_totallight <- ci_df %>%
  filter(variable == 'total_incoming_light') %>%
  select(-variable) 

fitted_totallight <- fitted_totallight %>%
  left_join(cf_totallight) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

pred_totallight <- pred_totallight %>%
  left_join(cf_totallight) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))


write.csv(pred_indivlight, file.path(fp_out, 'pred_indivlight.csv'), row.names = FALSE)
write.csv(pred_totallight, file.path(fp_out, 'pred_totallight.csv'), row.names = FALSE)
write.csv(fitted_indivlight, file.path(fp_out, 'fitted_indivlight.csv'), row.names = FALSE)
write.csv(fitted_totallight, file.path(fp_out, 'fitted_totallight.csv'), row.names = FALSE)

## Volume
ci_df <- read.csv('finalcsvs/volume_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

cf_volume <- read_csv('finalcsvs/volume_piecewise_cf_by_fg.csv') %>% 
  select(prod_model, fg, q50) %>% 
  rename(corr_factor = q50) %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))

fitted_totalvol <- ci_df %>%
  filter(variable == 'total_crown_volume_fitted') %>%
  select(-variable) 

fitted_totalvol <- fitted_totalvol %>%
  left_join(cf_volume) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

write.csv(fitted_totalvol, file.path(fp_out, 'fitted_totalvol.csv'), row.names = FALSE)

## Light captured
ci_df <- read.csv('finalcsvs/lightcaptured_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

cf_totallightcaptured <- read_csv('finalcsvs/lightcaptured_piecewise_cf_by_fg.csv') %>% 
  select(prod_model, fg, q50) %>% 
  rename(corr_factor = q50) %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

fitted_indivlightcaptured <- ci_df %>%
  filter(variable == 'captured_light_fitted') %>%
  select(-variable)

fitted_totallightcaptured <- ci_df %>%
  filter(variable == 'total_captured_light_fitted') %>%
  select(-variable) 

pred_indivlightcaptured <- ci_df %>%
  filter(variable == 'captured_light') %>%
  select(-variable)

pred_totallightcaptured <- ci_df %>%
  filter(variable == 'total_captured_light') %>%
  select(-variable) 

fitted_totallightcaptured <- fitted_totallightcaptured %>%
  left_join(cf_totallightcaptured) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

pred_totallightcaptured <- pred_totallightcaptured %>%
  left_join(cf_totallightcaptured) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))


write.csv(pred_indivlightcaptured, file.path(fp_out, 'pred_indivlightcaptured.csv'), row.names = FALSE)
write.csv(pred_totallightcaptured, file.path(fp_out, 'pred_totallightcaptured.csv'), row.names = FALSE)
write.csv(fitted_indivlightcaptured, file.path(fp_out, 'fitted_indivlightcaptured.csv'), row.names = FALSE)
write.csv(fitted_totallightcaptured, file.path(fp_out, 'fitted_totallightcaptured.csv'), row.names = FALSE)


# Diameter growth plots ---------------------------------------------------

# Observed individual diameter growth
obs_indivdiamgrowth_df <- map(alltreedat[-1],
                        function(x) with(x %>% filter(!recruit), cloudbin_across_years(dat_values = diam_growth_rate, dat_classes = dbh_corr, edges = binedgedata, n_census = 1)))
obs_indivdiamgrowth_df <- cbind(fg = 'all', year = rep(years, each = nrow(binedgedata)), do.call(rbind, obs_indivdiamgrowth_df))

obs_indivdiamgrowth_fg <- replicate(length(fg_names), list())

for (i in 1:length(fg_names)) {
  for (j in 1:length(years)) {
    obs_indivdiamgrowth_fg[[i]][[j]] <-
      cbind(fg = fg_names[i], year = years[j],
            cloudbin_across_years(dat_values = fgdat[[i]][[j+1]] %>% filter(!recruit) %>% pull(diam_growth_rate),
                                  dat_classes = fgdat[[i]][[j+1]] %>% filter(!recruit) %>% pull(dbh_corr),
                                  edges = binedgedata,
                                  n_census = 1))
  }
}

obs_indivdiamgrowth_fg <- do.call(rbind, map(obs_indivdiamgrowth_fg, function(x) do.call(rbind, x)))
obs_indivdiamgrowth_df <- rbind(obs_indivdiamgrowth_df, obs_indivdiamgrowth_fg)

# Fitted
fittedvals <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/diamgrowth_piecewise_ci_by_fg.csv'), stringsAsFactors = FALSE) %>%
  filter(variable == 'diameter_growth_fitted') %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))

# Added 09 Jan 2020: Correct bad numbers of individuals in observed diameter growth by replacing with values from density
obs_indivdiamgrowth_df$mean_n_individuals <- joined_counts$bin_count

write.csv(obs_indivdiamgrowth_df, file = file.path(fp_out, 'obs_indivdiamgrowth.csv'), row.names = FALSE)
write.csv(fittedvals, file = file.path(fp_out, 'fitted_indivdiamgrowth.csv'), row.names = FALSE)
