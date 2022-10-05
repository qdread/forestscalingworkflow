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

# crown volume and leaf area bins 
totalvolbins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = crownvolume, edges = binedgedata))
totalleafareabins_all <- with(alltree_light_95, logbin_setedges(x = dbh_corr, y = leaf_area, edges = binedgedata))

totalvolbins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$crownvolume, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totalvolbins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)
totalleafareabins_fg <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), 'unclassified')) %>%
  group_by(fg) %>%
  do(logbin_setedges(x = .$dbh_corr, y = .$leaf_area, edges = binedgedata)) %>%
  ungroup %>%
  rbind(data.frame(fg = 'all', totalleafareabins_all)) %>%
  mutate(bin_value = bin_value / area_core, year = 1995)

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
# total volume scaling

library(dplyr)

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

## Leaf area
ci_df <- read.csv('finalcsvs/leafarea_piecewise_ci_by_fg.csv', stringsAsFactors = FALSE)
area_core <- 42.84

ci_df$fg[ci_df$fg == 'alltree'] <- 'all'

cf_leafarea <- read_csv('finalcsvs/leafarea_piecewise_cf_by_fg.csv') %>% 
  select(prod_model, fg, q50) %>% 
  rename(corr_factor = q50) %>%
  mutate(fg = if_else(fg == 'alltree', 'all', fg))

fitted_totalleafarea <- ci_df %>%
  filter(variable == 'total_leaf_area_fitted') %>%
  select(-variable) 

fitted_totalleafarea <- fitted_totalleafarea %>%
  left_join(cf_leafarea) %>%
  mutate_at(vars(starts_with('q')), ~ . * (corr_factor/area_core))

write.csv(fitted_totalleafarea, file.path(fp_out, 'fitted_totalleafarea.csv'), row.names = FALSE)

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
