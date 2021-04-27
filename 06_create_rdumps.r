# ALL DATA DUMPS FOR (FINAL) FORESTLIGHT MODEL FITTING
# This script creates lists and writes them to Rdump text files that can be read into CmdStan for model fitting.
# QDR / ForestLight / 25 Oct 2019

# Modified 09 Jan 2020: remove recruits from the individual production input data.

# Needed data dumps:
# DBH and biomass production for density-biomass growth scalings
# DBH and diameter growth rate for density-diameter growth scalings
# Light per crown area and growth per crown area for growth/area-light/area scalings

# --------- #
# Load data #
# --------- #

# file path for writing rdumps (can be changed)
fpdump <- 'stanrdump'

# load raw data
load('data/rawdataobj_alternativecluster.r')

library(tidyverse)
library(rstan)

# read mortality data for mortality rdump
mort <- read_csv('data/data_forplotting/obs_mortalityindividuals.csv')

# -------------------------------- #
# Define function to create rdumps #
# -------------------------------- #

create_rdump <- function(dat, xvar, yvar, file_name = NULL, subsample = NULL, random_seed = NULL, x_range = c(1, 286)) {
  require(rstan)
  
  # Subset of data if specified
  if (!is.null(subsample) && nrow(dat) > subsample) {
	  if (!is.null(random_seed)) set.seed(random_seed)
    dat <- dat[sample(nrow(dat), subsample, replace = FALSE), ]
  }
  
  # Create list. Min and max x are taken from data, lower and upper limits are specified manually (the latter is a legacy from the Weibull fits)
  x <- dat[, xvar]
  y <- dat[, yvar]
  xdat <- list(N = length(x), x = x, y = y, x_min = min(x), x_max = max(x), LL = x_range[1], UL = x_range[2])
  
  if (!is.null(file_name)) {
    with(xdat, stan_rdump(names(xdat), file = file_name))
  } else {
    return(xdat)
  }
}

# -------------------------------- #
# Create rdumps and write to files #
# -------------------------------- #

# Do 1995 in all cases. For all cases, do all trees together and one for each FG
fgs <- c('fg1','fg2','fg3','fg4','fg5','unclassified')

### biomass growth ~ diameter scalings

create_rdump(alltreedat[[3]] %>% filter(!recruit), 'dbh_corr', 'production', file_name = file.path(fpdump, 'dump_production_alltree_1995.r'))

iwalk(fgs, ~ create_rdump(fgdat[[.y]][[3]] %>% filter(!recruit), 'dbh_corr' , 'production', file_name = file.path(fpdump, paste0('dump_production_', .x, '_1995.r'))))

### diameter growth ~ diameter scalings
create_rdump(alltreedat[[3]] %>% filter(!recruit), 'dbh_corr', 'diam_growth_rate', file_name = file.path(fpdump, 'dump_diamgrowthscaling_alltree_1995.r'))

iwalk(fgs, ~ create_rdump(fgdat[[.y]][[3]] %>% filter(!recruit), 'dbh_corr' , 'diam_growth_rate', file_name = file.path(fpdump, paste0('dump_diamgrowthscaling_', .x, '_1995.r'))))

### light per area ~ growth per area scalings

dat95 <- alltree_light_95 %>%
  select(dbh_corr, production, light_received, light_captured, crownarea, crownvolume, fg, recruit) %>%
  mutate(production_area = production/crownarea, light_area = light_received/crownarea) %>%
  mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg)))
  
create_rdump(dat95 %>% filter(!recruit), 'light_area', 'production_area', file_name = file.path(fpdump, 'dump_light_alltree_1995.r'))

dat95 %>%
  filter(!recruit) %>%
	group_by(fg) %>%
	group_walk(~ create_rdump(as.data.frame(.), 'light_area', 'production_area', file_name = file.path(fpdump, paste0('dump_light_', .y, '_1995.r'))))
	
### individual incoming light ~ diameter scalings
create_rdump(dat95, 'dbh_corr', 'light_received', file_name = file.path(fpdump, 'dump_rawlightscaling_alltree_1995.r'))

dat95 %>%
	group_by(fg) %>%
	group_walk(~ create_rdump(as.data.frame(.), 'dbh_corr', 'light_received', file_name = file.path(fpdump, paste0('dump_rawlightscaling_', .y, '_1995.r'))))

### individual captured light ~ diameter scalings (light received * fraction light captured; Added 26 Apr 2021)
create_rdump(dat95, 'dbh_corr', 'light_captured', file_name = file.path(fpdump, 'dump_lightcapturedscaling_alltree_1995.r'))

dat95 %>%
  group_by(fg) %>%
  group_walk(~ create_rdump(as.data.frame(.), 'dbh_corr', 'light_captured', file_name = file.path(fpdump, paste0('dump_lightcapturedscaling_', .y, '_1995.r'))))

### crown volume ~ diameter scalings
create_rdump(dat95, 'dbh_corr', 'crownvolume', file_name = file.path(fpdump, 'dump_volumescaling_alltree_1995.r'))

dat95 %>%
	group_by(fg) %>%
	group_walk(~ create_rdump(as.data.frame(.), 'dbh_corr', 'crownvolume', file_name = file.path(fpdump, paste0('dump_volumescaling_', .y, '_1995.r'))))

### leaf area ~ diameter scalings
create_rdump(dat95, 'dbh_corr', 'leaf_area', file_name = file.path(fpdump, 'dump_leafareascaling_alltree_1995.r'))

dat95 %>%
  group_by(fg) %>%
  group_walk(~ create_rdump(as.data.frame(.), 'dbh_corr', 'leaf_area', file_name = file.path(fpdump, paste0('dump_leafareascaling_', .y, '_1995.r'))))

### mortality
mort_data <- mort %>%
  filter(!fg %in% 'unclassified') %>%
  mutate(died = alive == 0) %>%
  select(fg, died, light_received_byarea)

mort_data_dump <- with(mort_data, list(N = nrow(mort_data), M = 5, x = light_received_byarea, y = as.numeric(died), fg = as.numeric(factor(fg))))

with(mort_data_dump, stan_rdump(list = names(mort_data_dump), file = file.path(fpdump, 'mortalitydump.r')))

# ------------------------------------------------------------------------------ #
# Create CSVs of minima, maxima, number of individuals, and normalization totals #
# ------------------------------------------------------------------------------ #

# This info is used to make some corrections in the output after the models have been run, for plotting purposes.

years <- c(1985, 1990, 1995, 2000, 2005, 2010)

# Minima, maxima, and number of individuals 
valall <- map2_dfr(alltreedat, years,
				   ~ data.frame(fg = 'alltree', year = .y, xmin = min(.x$dbh_corr), n = nrow(.x), n_not_recruit = sum(!.x$recruit)))

valfg <- map2_dfr(fgdat, fgs, function(dat, fg_name) {
	map2_dfr(dat, years,
			 ~ data.frame(fg = fg_name, year = .y, xmin = min(.x$dbh_corr), n = nrow(.x), n_not_recruit = sum(!.x$recruit)))
})

min_n <- rbind(valall, valfg)
write_csv(min_n, file.path(fpdump, 'min_n.csv'))

# Minima, maxima, and number of individuals for 1995 only for trees with light measured
valall <- data.frame(fg = 'alltree', year = 1995, xmin = with(alltree_light_95, min(dbh_corr)), n = nrow(alltree_light_95), n_not_recruit = sum(!alltree_light_95$recruit))
valfg <- alltree_light_95 %>%
  group_by(fg) %>%
  summarize(xmin = min(dbh_corr),
            n = n(),
            n_not_recruit = sum(!recruit))

min_n <- data.frame(fg = c('alltree', 'fg1','fg2','fg3','fg4','fg5','unclassified'),
                    year = 1995,
                    xmin = c(valall$xmin, valfg$xmin),
                    n = c(valall$n, valfg$n))

write_csv(min_n, file.path(fpdump, 'min_n_lighttrees.csv'))

