# Script 3: Create log-bins for all variables.

# Modified on 01 April 2020: change nomenclature
# Modified on 20 February 2020: For ratios, use minimum number of individuals across the 2 groups being compared, rather than total
# Modified on 13 January 2020: Use imputed production (1 segment) to create total production bins -- note this must now be run AFTER models are all fit.
# Modified on 13 January 2020: For individual production, diameter growth, and production per area, remove the recruits.
# Modified on 06 November 2019: prioritize single year binning, separating multiple year binning elsewhere
# Modified on 05 November 2019: Use packaged functions
# Script split off on 30 October 2019.
# Crown volume added on 21 March 2019.

# Procedure: use all 5 FGs (but not unclassified) together to get the bin edges
# Then apply those bin edges to each FG in isolation.
# That way each FG are given the same bin edges.

library(tidyverse)
library(forestscaling)

load('data/rawdataobj_withimputedproduction.RData')
group_names <- c('all','all_classified','fg1','fg2','fg3','fg4','fg5','unclassified')
fg_names <- c('fg1','fg2','fg3','fg4','fg5','unclassified')
years <- seq(1990, 2010, 5)


# Set number of bins       
numbins <- 20

# Make a version of alltreedat without the unclassified trees
alltreedat_classified <- map(alltreedat, ~ filter(., !is.na(fg)))

# Bin classified trees. (log binning of density)
allyeardbh_classified <- map(alltreedat_classified[-1], ~ pull(., dbh_corr)) %>% unlist
dbhbin_allclassified <- logbin(x = allyeardbh_classified, y = NULL, n = numbins)

# Bin all trees including unclassified
allyeardbh <- map(alltreedat[-1], ~ pull(., dbh_corr)) %>% unlist
dbhbin_all <- logbin(x = allyeardbh, y = NULL, n = numbins)

allyeardbh_fg <- fgdat %>% map(~ map(., ~ pull(., dbh_corr)) %>% unlist)

# The dbhbin_ objects can be used as the edges argument in logbin_setedges()

# Make versions of input data without the new recruits
alltreedat_norecruits <- map(alltreedat, ~ filter(., !recruit))
alltreedat_classified_norecruits <- map(alltreedat_classified, ~ filter(., !recruit))
fgdat_norecruits <- map(fgdat, ~ map(., ~ filter(., !recruit)))

allyeardbh_norecruits <- map(alltreedat_norecruits[-1], ~ pull(., dbh_corr)) %>% unlist
allyeardbh_classified_norecruits <- map(alltreedat_classified_norecruits[-1], ~ pull(., dbh_corr)) %>% unlist
allyeardbh_fg_norecruits <- fgdat_norecruits %>% map(~ map(., ~ pull(., dbh_corr)) %>% unlist)

# Bin density and individual production by year using the above generated bin edges.

dbhbin_all_byyear <- alltreedat[-1] %>% map(~ logbin_setedges(x = .$dbh_corr, y = NULL, edges = dbhbin_all))
dbhbin_allclassified_byyear <- alltreedat_classified[-1] %>% map(~ logbin_setedges(x = .$dbh_corr, y = NULL, edges = dbhbin_allclassified))

dbhbin_fg_byyear <- fgdat %>%
  map(~ map(.[-1], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_allclassified)))

dbhbin_all_5census <- bin_across_years(dbhbin_all_byyear)
dbhbin_allclassified_5census <- bin_across_years(dbhbin_allclassified_byyear)
dbhbin_fg_5census <- lapply(dbhbin_fg_byyear, bin_across_years)

# Combine single year bin into single df
densitybin_byyear <- rbind(
  data.frame(fg = 'all', map2_dfr(dbhbin_all_byyear, years, ~ data.frame(year = .y, .x))),
  map2_dfr(dbhbin_fg_byyear, fg_names, ~ data.frame(fg = .y, map2_dfr(.x, years, ~ data.frame(year = .y, .x))))
)

# Combine multiple year bin into single df
densitybin_5census <- cbind(fg = rep(group_names, each = numbins), 
                            rbind(dbhbin_all_5census, dbhbin_allclassified_5census, do.call('rbind', dbhbin_fg_5census)))

# Individual production 
# Take the mean and 2.5, 50 (median), 97.5 quantiles within each bin.
# Do it across all years, and for separate years.
# Exclude new recruits (added 13 Jan 2020)
prodbin_all_byyear <- alltreedat_norecruits[-1] %>% map(~ cloudbin_across_years(dat_values = .$production, dat_classes = .$dbh_corr, edges = dbhbin_allclassified, n_census = 1))
prodbin_fg_byyear <- alltreedat_norecruits[-1] %>%
  map(function(dat) {
    dat %>%
      group_by(fg) %>%
      do(cloudbin_across_years(dat_values = .$production, dat_classes = .$dbh_corr, edges = dbhbin_allclassified, n_census = 1))
  })

allyearprod <- map(alltreedat_norecruits[-1], ~ pull(., production)) %>% unlist
allyearprod_classified <- map(alltreedat_classified_norecruits[-1], ~ pull(., production)) %>% unlist
allyearprod_fg <- fgdat_norecruits %>% map(~ map(., ~ pull(., production)) %>% unlist)

# Combine single year bin into single df.
indivproductionbin_byyear <- bind_rows(
  data.frame(fg = 'all', map2_dfr(prodbin_all_byyear, years, ~ data.frame(year = .y, .x))),
  map2_dfr(prodbin_fg_byyear, years, ~ data.frame(year = .y, .x) %>% mutate(fg = if_else(is.na(fg), 'unclassified', paste0('fg', fg))))
) %>%
  select(-mean_n_individuals)

# Combine multiple year bin into single df.
prodbin_all_5census <- cloudbin_across_years(dat_values = allyearprod, dat_classes = allyeardbh_norecruits, edges = dbhbin_all)
prodbin_allclassified_5census <- cloudbin_across_years(dat_values = allyearprod_classified, dat_classes = allyeardbh_classified_norecruits, edges = dbhbin_allclassified)

prodbin_fg_5census <- map2(allyearprod_fg, allyeardbh_fg_norecruits, ~ cloudbin_across_years(dat_values = .x, dat_classes = .y, edges = dbhbin_allclassified))

indivproductionbin_5census <- cbind(fg = rep(group_names, each = numbins),
                                    rbind(prodbin_all_5census, prodbin_allclassified_5census, do.call('rbind', prodbin_fg_5census)))

# Total production
# Do the binning for each year separately, as for density, then find min, max, and median.

totalprodbin_alltree_byyear <- alltreedat[-1] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$production_imputed1, edges = dbhbin_all))
totalprodbin_allclassified_byyear <- alltreedat_classified[-1] %>% map(~ logbin_setedges(x = .$dbh_corr, y = .$production_imputed1, edges = dbhbin_allclassified))

totalprodbin_fg_byyear <- fgdat %>%
  map(~ map(.[-1], function(z) logbin_setedges(x = z$dbh_corr, y = z$production_imputed1, edges = dbhbin_allclassified)))

totalprodbin_all_5census <- bin_across_years(totalprodbin_alltree_byyear)
totalprodbin_allclassified_5census <- bin_across_years(totalprodbin_allclassified_byyear)
totalprodbin_fg_5census <- lapply(totalprodbin_fg_byyear, bin_across_years)

# Combine single year bin into single df
totalproductionbin_byyear <- rbind(
  data.frame(fg = 'all', map2_dfr(totalprodbin_alltree_byyear, years, ~ data.frame(year = .y, .x))),
  map2_dfr(totalprodbin_fg_byyear, fg_names, ~ data.frame(fg = .y, map2_dfr(.x, years, ~ data.frame(year = .y, .x))))
)

# Combine multiple year bin into single df
totalproductionbin_5census <- cbind(fg = rep(group_names, each = numbins), 
                                    rbind(totalprodbin_all_5census, totalprodbin_allclassified_5census, do.call('rbind', totalprodbin_fg_5census)))

# Figure 5 ratio binning
# Previously we had shade:gap ratio of both production and density by light received bin, as well as average shade score for each light bin.
# This was repeated for diameter bins.
# Here do the same 6 binning types, but repeat twice, once for ratio fg2:fg4 (longlived pioneer to shortlived breeder)
# and once for ratio fg1:fg3 (fast:slow)

totalprodbin_byyear_bydiam <- fgdat %>%
  map(~ map(.[-1], function(z) logbin_setedges(x = z$dbh_corr, y = z$production_imputed1, edges = dbhbin_allclassified)))

densitybin_byyear_bydiam <- fgdat %>%
  map(~ map(.[-1], function(z) logbin_setedges(x = z$dbh_corr, y = NULL, edges = dbhbin_allclassified)))

# pioneer to breeder by diameter (fg2 to fg4)
breeder_stats_bydiam <- tibble(fg_a_prod = totalprodbin_byyear_bydiam[[2]],
                               fg_b_prod = totalprodbin_byyear_bydiam[[4]],
                               fg_a_dens = densitybin_byyear_bydiam[[2]],
                               fg_b_dens = densitybin_byyear_bydiam[[4]],
                               year = c(1990,1995,2000,2005,2010)) %>%
  pmap(function(fg_a_prod, fg_b_prod, fg_a_dens, fg_b_dens, year) 
    data.frame(bin = 1:numbins,
               year = year,
               breeder_production_ratio = fg_a_prod$bin_value / fg_b_prod$bin_value,
               breeder_density_ratio = fg_a_dens$bin_value / fg_b_dens$bin_value,
               n_individuals = pmin(fg_a_dens$bin_count, fg_b_dens$bin_count))) %>%
  bind_rows %>%
  mutate_if(is.double, ~ if_else(is.finite(.x), .x, as.numeric(NA)))

breeder_stats_bydiam_2census <- breeder_stats_bydiam %>% 
  filter(year %in% c(1990, 1995)) %>%
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(breeder_production_ratio),
            density_ratio_mean = mean(breeder_density_ratio),
            production_ratio_min = min(breeder_production_ratio),
            production_ratio_max = max(breeder_production_ratio),
            density_ratio_min = min(breeder_density_ratio),
            density_ratio_max = max(breeder_density_ratio),
            mean_n_individuals = min(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

breeder_stats_bydiam_byyear <- breeder_stats_bydiam %>%
  select(-bin) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

breederscore_bin_bydiam_2census <- binscore(dat = alltreedat[2:3], bindat = dbhbin_allclassified, score_column = 'X2', class_column = 'dbh_corr')
breederscore_bin_bydiam_byyear <- map2_dfr(alltreedat[2:6], years, ~ data.frame(year = .y, cloudbin_across_years(dat_values = .x$X2[!is.na(.x$X2)], dat_classes = .x$dbh_corr[!is.na(.x$X2)], edges = dbhbin_allclassified, mean = 'arithmetic', n_census = 1))) %>%
  select(-mean_n_individuals)

breeder_stats_bydiam_5census <- breeder_stats_bydiam %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(breeder_production_ratio),
            density_ratio_mean = mean(breeder_density_ratio),
            production_ratio_min = min(breeder_production_ratio),
            production_ratio_max = max(breeder_production_ratio),
            density_ratio_min = min(breeder_density_ratio),
            density_ratio_max = max(breeder_density_ratio),
            mean_n_individuals = min(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 


# Fast to slow by diameter (fg1 to fg3)
fastslow_stats_bydiam <- tibble(fg_a_prod = totalprodbin_byyear_bydiam[[1]],
                                fg_b_prod = totalprodbin_byyear_bydiam[[3]],
                                fg_a_dens = densitybin_byyear_bydiam[[1]],
                                fg_b_dens = densitybin_byyear_bydiam[[3]],
                                year = c(1990,1995,2000,2005,2010)) %>%
  pmap(function(fg_a_prod, fg_b_prod, fg_a_dens, fg_b_dens, year) 
    data.frame(bin = 1:numbins,
               year = year,
               fastslow_production_ratio = fg_a_prod$bin_value / fg_b_prod$bin_value,
               fastslow_density_ratio = fg_a_dens$bin_value / fg_b_dens$bin_value,
               n_individuals = pmin(fg_a_dens$bin_count, fg_b_dens$bin_count))) %>%
  bind_rows %>%
  mutate_if(is.double, ~ if_else(is.finite(.x), .x, as.numeric(NA)))

fastslow_stats_bydiam_byyear <- fastslow_stats_bydiam %>%
  select(-bin) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

fastslow_stats_bydiam_2census <- fastslow_stats_bydiam %>% 
  filter(year %in% c(1990, 1995)) %>%
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(fastslow_production_ratio),
            density_ratio_mean = mean(fastslow_density_ratio),
            production_ratio_min = min(fastslow_production_ratio),
            production_ratio_max = max(fastslow_production_ratio),
            density_ratio_min = min(fastslow_density_ratio),
            density_ratio_max = max(fastslow_density_ratio),
            mean_n_individuals = min(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

fastslowscore_bin_bydiam_2census <- binscore(dat = alltreedat[2:3], bindat = dbhbin_allclassified, score_column = 'X1', class_column = 'dbh_corr')
fastslowscore_bin_bydiam_byyear <- map2_dfr(alltreedat[2:6], years, ~ data.frame(year = .y, cloudbin_across_years(dat_values = .x$X1[!is.na(.x$X1)], dat_classes = .x$dbh_corr[!is.na(.x$X1)], edges = dbhbin_allclassified, mean = 'arithmetic', n_census = 1))) %>%
  select(-mean_n_individuals)


fastslow_stats_bydiam_5census <- fastslow_stats_bydiam %>% 
  group_by(bin) %>%
  summarize(production_ratio_mean = mean(fastslow_production_ratio),
            density_ratio_mean = mean(fastslow_density_ratio),
            production_ratio_min = min(fastslow_production_ratio),
            production_ratio_max = max(fastslow_production_ratio),
            density_ratio_min = min(fastslow_density_ratio),
            density_ratio_max = max(fastslow_density_ratio),
            mean_n_individuals = min(n_individuals)) %>%
  cbind(densitybin_byyear_bydiam[[1]][[1]][,c('bin_midpoint', 'bin_min', 'bin_max')]) 

# Export binned data

fpdata <- 'data/data_binned'

# Multiple year bins
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'breeder_stats_bydiam_2census',  'breederscore_bin_bydiam_2census',  'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census','fastslow_stats_bydiam_5census','breeder_stats_bydiam_5census')

for (i in file_names) {
  write.csv(get(i), file=file.path(fpdata, paste0(i,'.csv')), row.names = FALSE)
}

save(list = file_names, file = file.path(fpdata, 'bin_object_multipleyear.RData'))

# Single year bins
file_names <- c('densitybin_byyear', 'indivproductionbin_byyear', 'totalproductionbin_byyear', 'breeder_stats_bydiam_byyear', 'breederscore_bin_bydiam_byyear', 'fastslow_stats_bydiam_byyear', 'fastslowscore_bin_bydiam_byyear') 

for (i in file_names) {
  write.csv(get(i), file=file.path(fpdata, paste0(i,'.csv')), row.names = FALSE)
}

save(list = file_names, file = file.path(fpdata, 'bin_object_singleyear.RData'))

