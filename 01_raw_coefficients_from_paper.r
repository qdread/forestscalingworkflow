# Script to extract table S2 from Martinez-Cano et al. 2019 and convert to a clean CSV (allometric coefficients for tree height and crown area)
# Tree height in m, crown area in m^2, from dbh in cm.
# QDR / ForestLight / 11 Nov 2019

library(tabulizer)
library(tidyverse)

martinezS2raw <- extract_tables('data/bg-16-847-2019-supplement.pdf', pages = 4:7) # Table S2 only

map(martinezS2raw, head) # All have the same two header rows

# Lop off headers and join to single data frame
martinezS2 <- map_dfr(martinezS2raw, ~ .[-(1:2),] %>% as.data.frame %>% setNames(c('species', 'binomial', 'height_pars', 'area_pars')))

# Fix the one species name that ran onto two rows
rowidx <- grep('Protium tenuifolium', martinezS2$binomial)

martinezS2$binomial[rowidx + 1] <- 'Protium tenuifolium subsp. sessiliflorum'
martinezS2 <- martinezS2[-rowidx, ]


# Separate parameter estimates and credible interval bounds into their own columns
height_par_names <- paste('height', c('a', 'a_min', 'a_max', 'b', 'b_min', 'b_max', 'k', 'k_min', 'k_max'), sep = '_')
area_par_names <- paste('area', c('a', 'a_min', 'a_max', 'b', 'b_min', 'b_max'), sep = '_')

strtopardf <- function(pars, par_names) {
  extr <- str_extract_all(pars, '[0-9.-]+')
  do.call(rbind, extr) %>%
    apply(2, as.numeric) %>%
    as.data.frame %>%
    setNames(par_names)
}

height_pars_extracted <- strtopardf(martinezS2$height_pars, height_par_names)
area_pars_extracted <- strtopardf(martinezS2$area_pars, area_par_names)

martinezS2 <- martinezS2 %>%
  select(species, binomial) %>%
  cbind(height_pars_extracted, area_pars_extracted)

write_csv(martinezS2, 'data/allometry_martinez-cano_table_s2.csv')
