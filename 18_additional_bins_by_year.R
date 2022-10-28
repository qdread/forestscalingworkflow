# Create binned data for 5 censuses to be used to show range of min-max binned values by year for each diameter class bin in each FG
# QDR 2022-10-27

library(tidyverse)
library(forestscaling)

load('data/rawdataobj_alternativecluster.r')

# Load 1995 bin data (to make sure we're using the same bin breaks as other figs)
load('data/data_binned/bin_object_singleyear.RData')

# Create binned richness data

# Use existing bin bounds
bin_bounds <- fastslow_stats_bydiam_byyear[1:20, c('bin_midpoint', 'bin_min', 'bin_max')]

# Set up bin data structure.
years <- seq(1990, 2010, by=5)

bin_fg_year <- expand_grid(year = years, fg = c(paste0('fg', 1:5), 'unclassified'), bin = 1:20) %>%
  left_join(bin_bounds %>% mutate(bin = 1:20))

# Concatenate data from the five relevant censuses 1990-2010
dat <- map2_dfr(alltreedat[-1], years, function(x,y) {
  x %>%
    select(sp, fg, dbh_corr, DFstatus, production) %>%
    mutate(year = y)
}) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))

### RICHNESS BINS (ALSO INCLUDES RICHNESS VS. ABUNDANCE)
richness_bin_fg_year <- bin_fg_year %>%
  cbind(pmap_dfr(bin_fg_year, function(year, fg, bin_min, bin_max, ...) {
    dat_subset <- dat[dat$fg %in% fg & dat$year %in% year & dat$dbh_corr >= bin_min & dat$dbh_corr < bin_max, ]
    sp_ids <- as.character(dat_subset$sp)
    data.frame(richness = length(unique(sp_ids)),
               n_individuals = length(sp_ids),
               production = sum(dat_subset$production))
  })) %>%
  mutate(richness_by_bin_width = richness / (bin_max - bin_min),
         abundance_by_bin_width = n_individuals / (bin_max - bin_min))

### RICHNESS & ABUNDANCE RATIO BINS
# Use the previous result and divide fast:slow (1/3) and pioneer:slow (2/3)

bin_ratio_year <- richness_bin_fg_year %>%
  filter(fg %in% c('fg1','fg2','fg3')) %>%
  pivot_wider(id_cols = c(year, bin_midpoint, bin_min, bin_max), names_from = fg, values_from = c(richness, production, n_individuals)) %>%
  mutate(richness_ratio_fastslow = richness_fg1/richness_fg3,
         richness_ratio_pioneerslow = richness_fg2/richness_fg3,
         abundance_ratio_fastslow = n_individuals_fg1/n_individuals_fg3,
         abundance_ratio_pioneerslow = n_individuals_fg2/n_individuals_fg3,
         production_ratio_fastslow = production_fg1/production_fg3,
         production_ratio_pioneerslow = production_fg2/production_fg3,
         min_n_individuals_fastslow = pmin(n_individuals_fg1, n_individuals_fg3),
         min_n_individuals_pioneerslow = pmin(n_individuals_fg2, n_individuals_fg3)) %>%
  select(year, bin_midpoint, bin_min, bin_max, contains('fast'), contains('pioneer'))

### MORTALITY BINS
# This will require finding the mortality rate of each year.
# Make a column that says which year was the last time that treeID was in the census.
# The "died" column will indicate if the tree died in the interval following the census.
dat_mort <- map2_dfr(alltreedat, c(1985, years), function(x,y) {
  x %>%
    select(treeID, sp, fg, dbh_corr) %>%
    mutate(year = y)
}) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))) %>%
  group_by(treeID) %>%
  mutate(last_year = max(year)) %>% 
  ungroup %>%
  mutate(died = year == last_year)

bin_fg_year <- expand_grid(year = seq(1985, 2005, by = 5), fg = c(paste0('fg', 1:5), 'unclassified'), bin = 1:20) %>%
  left_join(bin_bounds %>% mutate(bin = 1:20))

bin_mort_year <- bin_fg_year %>%
  mutate(mortality = pmap_dbl(bin_fg_year, function(year, fg, bin_min, bin_max, ...) {
    dat_mort_subset <- dat_mort[dat_mort$fg %in% fg & dat_mort$year %in% year & dat_mort$dbh_corr >= bin_min & dat_mort$dbh_corr < bin_max, ]
    sum(dat_mort_subset$died)/nrow(dat_mort_subset)
  })) %>%
  mutate(year = year + 5)

# Join mortality bin to the other by-fg bin.
additional_bins_fg_year <- richness_bin_fg_year %>% left_join(bin_mort_year)

write_csv(additional_bins_fg_year %>% rename(abundance = n_individuals), 'data/data_binned/additional_bins_fg_year.csv')
write_csv(bin_ratio_year, 'data/data_binned/additional_bins_ratio_year.csv')
