#--------Plots all main and supplemental figures in Grady et al's "Life history scaling in a tropical forest" ---------

########### ============== ##############
### WHICH MODEL FITS TO USE IN PLOTS ####
########### ============== ##############

# change these if needed.
DENS = 3
PROD = 1

# size of plot = 42.84 hectares

library(broom)
library(egg)
library(scales)
library(RColorBrewer)
library(gtable)
library(grid)
library(reshape2)
library(hexbin)
library(Hmisc, pos = 100)
library(rstan)
library(lme4)
library(sjstats)
library(rstanarm)
library(forestscaling) # Packaged functions and ggplot2 themes here! 
library(tidyverse)
library(ggnewscale)

# Define color schemes and labels
guild_fills_all = c(fg1 = "#BFE046", fg2 =  "#267038", fg3 = "#27408b", fg4 = "#87Cefa", fg5 = "gray93", all = "black")
guild_fills_fg = c(fg1 = "#BFE046", fg2 =  "#267038", fg3 = "#27408b", fg4 = "#87Cefa", fg5 = "gray93")

guild_colors_all <- c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray")
guild_colors_fg <- c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray")

slow_fast_fill <- c("#27408b", "#BFE046" ) 
slow_tall_fill <- c("#27408b", "#267038" ) #'tall' is long-lived pioneer
short_tall_fill <- c("#87Cefa", "#267038")  

scale_slow_tall <- scale_fill_gradientn(colors = slow_tall_fill, trans = 'log')
scale_slow_fast <- scale_fill_gradientn(colors = slow_fast_fill, trans = 'log')

scale_short_tall  <- scale_fill_gradientn(colors = short_tall_fill) #'short' is short-lived recruiter
scale_slow_fast2 <- scale_fill_gradientn(colors = slow_fast_fill) #no log for PCA values

fg_labels <- c("fg1", "fg2","fg3", "fg4", "fg5")

guild_lookup <- data.frame(fg = c('fg1','fg2','fg3','fg4','fg5','all','unclassified'), fg_name = c('Fast','LL Pioneer','Slow','Short','Medium','All','Unclassified'))
guild_lookup
year_to_plot = 1995
geom_size <- 4

grob_fast <- grobTree(textGrob("Fast", x = 0.04, y = 0.95,  hjust = 0, gp = gpar(col = "#BFE046", fontsize = 15, fontface = "italic"))) 
grob_tall <- grobTree(textGrob("LL Pioneer", x = 0.04, y = 0.88,  hjust = 0, gp = gpar(col = "#267038", fontsize = 15, fontface = "italic"))) 
grob_medium <- grobTree(textGrob("Medium", x = 0.04, y = 0.81,  hjust = 0, gp = gpar(col = "gray70", fontsize = 15, fontface = "italic"))) 
grob_slow <- grobTree(textGrob("Slow", x = 0.04, y = 0.74,  hjust = 0, gp = gpar(col = "#27408b", fontsize = 15, fontface = "italic"))) 
grob_short <- grobTree(textGrob("Short", x = 0.04, y = 0.67,  hjust = 0, gp = gpar(col = "#87Cefa", fontsize = 15, fontface = "italic"))) 
grob_all <- grobTree(textGrob("All", x = 0.04, y = 0.60,  hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "italic"))) 


# To add back the legend
theme_plant2 <- forestscaling::theme_plant() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

#adjusting theme
theme_plant <- forestscaling::theme_plant() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position = "right", legend.text = element_text(size = 14 ), legend.key = element_blank())

#---------------------------- Data ---------------------------- 
fp <-'data/data_forplotting'

for (i in dir(fp, pattern = 'obs_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}


for (i in dir(fp, pattern = 'pred_|fitted_')) {
  n <- gsub('.csv','',i)
  assign(n, read.csv(file.path(fp, i), stringsAsFactors = FALSE))
}
#### Extract model output to get the fitted values, slopes, etc.
load('data/data_piecewisefits/fits_bylight_forratio.RData')
load('data/data_binned/bin_object_singleyear.RData')
load('data/data_forplotting/light_scaling_plotting_data.RData')
load('data/rawdataobj1995.RData')
binned_data <- read_csv("data/data_binned/additional_bins_fg_year.csv")
ratios <- read_csv('data/data_binned/additional_bins_ratio_year.csv')
fgbci <- read.table('data/fgroups_dynamics_new.txt', stringsAsFactors = FALSE)
params <- read.csv('data/data_piecewisefits/piecewise_paramci_by_fg.csv', stringsAsFactors = FALSE)
growth_slopes <- read.csv('data/data_piecewisefits/piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE)
light_slopes <- read.csv('data/data_piecewisefits/light_piecewise_fitted_slopes_by_fg.csv', stringsAsFactors = FALSE)
totallightbins_fg <- read.csv(file.path(fp, 'obs_totallight.csv'), stringsAsFactors = FALSE)
unscaledlightbydbhcloudbin_fg <- read.csv(file.path(fp, 'unscaledlightbydbhcloudbin_fg.csv'), stringsAsFactors = FALSE)

# source the extra extraction functions that aren't in the package
source('R_functions/model_output_extraction_functions.r')
source('R_functions/extra_plot_functions.R') # add new plotting themes


###################################################################################
################# Fig 2: Abundance, Richness, Production   #####################
###################################################################################

#----------------------------------------------------------------------------------
#--------------------------   Fig 2A: Abundance          ----------------------------------
#----------------------------------------------------------------------------------

#--error bars
fg_mid = c("fg", "bin_midpoint")
dens_range <- obs_dens %>%
  filter(fg != "unclassified",bin_count >= 20) %>%
  group_by(fg, bin_midpoint) %>%
  summarize(
    min = min(bin_value),
    max = max(bin_value)
  )
dens_range

(p <- plot_dens(year_to_plot = 1995,
                fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                model_fit = DENS,
                dodge_width = 0.03,
                x_limits = c(.9, 160),
                y_limits = c(0.007, 20000),
                x_breaks = c(1, 10, 100),
                y_labels = c(0.001, 0.1, 10,1000),
                y_breaks = c(0.001, 0.1,  10, 1000)))
p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


#----------------------------------------------------------------------------------
#--------------------------   Fig S3: Abundance across census years ---------------
#----------------------------------------------------------------------------------
(p <- plot_dens2(year_to_plot = 1995,
                fg_names = c('fg1','fg2','fg3','fg4','fg5','all'),
                model_fit = DENS,
                dodge_width = 0.03,
                x_limits = c(.9, 160),
                y_limits = c(0.007, 20000),
                x_breaks = c(1, 10, 100),
                y_labels = c(0.001, 0.1, 10,1000),
                y_breaks = c(0.001, 0.1,  10, 1000)))

p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

#----------------------------------------------------------------------------------
#--------------------------   Fig. 2C Richness  ------------------------------------
#----------------------------------------------------------------------------------

# filter fitted sizes by binned sizes
max_dbh_fg <- obs_richnessbydiameter  %>%
  filter(n_individuals >= 20) %>%
  group_by(fg) %>%
  summarize(max_dbh = max(bin_midpoint) + 2 )

min_dbh_fg <-obs_richnessbydiameter  %>%
  filter(n_individuals >= 20) %>%
  group_by(fg) %>%
  summarize(min_dbh = min(bin_midpoint) - 0.1)

fitted_richnessbydiameter_filtered <- fitted_richnessbydiameter %>%
  left_join(max_dbh_fg) %>%
  left_join(min_dbh_fg) %>%
  group_by(fg) %>%
  filter(dbh <= max_dbh, dbh >= min_dbh)

unique(fitted_richnessbydiameter_filtered$year)

fg_mid = c("fg", "bin_midpoint")

rich_range <- binned_data %>%
  filter(fg != "unclassified", abundance >= 20) %>%
  group_by(across(all_of(fg_mid))) %>%
  summarize(
    min = min(richness_by_bin_width),
    max = max(richness_by_bin_width)
  )


#---- Plot Richness ------
(p <- ggplot() + 
    theme_plant() +
    scale_x_log10(name = 'Stem Diameter (cm)', limit = c(.9, 160)) + 
    scale_y_log10(labels = signif, limits = c(0.007, 20000), position = "left",
                  name = expression(paste("Richness (cm"^-1, ")"))) +
    scale_fill_manual(values = guild_fills_all) +
    scale_color_manual(values = guild_colors_all) +
    geom_ribbon(data = fitted_richnessbydiameter_filtered  %>% 
                  arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))),
                aes(x = dbh, ymin = q025, ymax = q975,
                    group = fg, fill = fg), alpha = 0.4) +
    geom_line(data = fitted_richnessbydiameter_filtered  %>% 
                arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))),
              aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_errorbar(data = rich_range %>%  arrange(desc(fg)), aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), width = 0) +
    geom_point(data = binned_data %>% 
                 arrange(desc(fg)) %>%
                 filter(!fg %in% 'unclassified', abundance >= 20) %>% 
                 filter(year == "1995") %>%
                 arrange(desc(fg)), 
               aes(x = bin_midpoint, y = richness_by_bin_width, fill = fg, color = fg),
               shape = 21, size = 4, color = "black") + theme_no_x() 
)

p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


#----------------------------------------------------------------------------------
#--------------------------   Fig. S3 Richness  across census years ---------------
#----------------------------------------------------------------------------------

(p <- ggplot() + 
   theme_plant() +
   scale_x_log10(name = 'Stem Diameter (cm)',
                 limit = c(.9, 160)) + 
   scale_y_log10(labels = signif,
                 limits = c(0.007, 20000),
                 position = "left",
                 name = expression(paste("Richness (cm"^-1, ")"))) +
   scale_fill_manual(values = guild_fills_all) +
   scale_color_manual(values = guild_colors_all) +
   geom_line(data = fitted_richnessbydiameter_filtered  %>% 
               arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))),
             aes(x = dbh, y = q50, group = fg, color = fg), size = 0.2) +
   geom_errorbar(data = rich_range %>%  arrange(desc(fg)), aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), size = .5, width = 0) +
   theme_no_x() 
)

p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


#----------------------------------------------------------------------------------
#--------------------------   Fig. 2E. Productivity  ------------------------------
#----------------------------------------------------------------------------------

obs_totalprod <- obs_totalprod %>%
  filter(bin_count >= 20)

fg_mid = c("fg", "bin_midpoint")

prod_range <- obs_totalprod %>%
  group_by(across(all_of(fg_mid))) %>%
  filter(fg != "unclassified", bin_count >= 20) %>%
  summarize(
    min = min(bin_value),
    max = max(bin_value)
  )

p <- plot_totalprod(year_to_plot = 1995,
                    fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                    model_fit_density = DENS, 
                    model_fit_production = PROD,
                    x_name = "Stem Diameter (cm)", 
                    x_limits = c(0.9,160),
                    y_limits = c(0.3, 200),
                    y_breaks = c(0.1, 1, 10, 100),
                    y_labels = c(0.1, 1, 10, 100),
                    dodge_width = 0.03,
                    preddat = fitted_totalprod)


p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


#----------------------------------------------------------------------------------
#-------------------  Fig. S#. Productivity across census years -------------------
#----------------------------------------------------------------------------------

p <- plot_totalprod2(year_to_plot = 1995,
                    fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                    model_fit_density = DENS, 
                    model_fit_production = PROD,
                    x_name = "Stem Diameter (cm)", 
                    x_limits = c(0.9,160),
                    y_limits = c(0.3, 200),
                    y_breaks = c(0.1, 1, 10, 100),
                    y_labels = c(0.1, 1, 10, 100),
                    dodge_width = 0.03,
                    preddat = fitted_totalprod)


p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


########################################################################
######################   Fig 2 B,D,E: Ratios     #######################
########################################################################



#----------------------------------------------------------------------------------
#--------------------------   Fig. 2B Relative Abundance  -------------------------
#----------------------------------------------------------------------------------

abun_ratio_pioneerslow_range <- ratios %>%
  group_by(bin_midpoint) %>%
  filter(min_n_individuals_pioneerslow >= 20) %>%
  summarize(
    min = min(abundance_ratio_pioneerslow),
    max = max(abundance_ratio_pioneerslow))

abun_ratio_fastslow_range <- ratios %>%
  group_by(bin_midpoint) %>%
  filter(min_n_individuals_fastslow >= 20) %>%
  summarize(
    min = min(abundance_ratio_fastslow),
    max = max(abundance_ratio_fastslow))

abun_data = ratios %>% filter(min_n_individuals_pioneerslow >= 20, min_n_individuals_fastslow >= 20, year == "1995")

abun_max = max(abun_data$bin_midpoint)

#------- Fig 2B-  
(abun_ratio  <- ggplot() + 
    theme_plant() + theme_no_x() + 
    scale_x_log10(limits = c(.9, 100), breaks = c(1, 10, 100), 
                  position = "bottom", 
                  expression(paste('Stem Diameter (cm)'))) + 
    scale_y_log10(labels = signif, 
                  breaks = c(0.01, 0.1, 1, 10), limits = c(0.03, 10),
                  position = "right", name = expression("Abundance Ratio")) + 
    stat_smooth(data = ratios %>% filter(min_n_individuals_pioneerslow >= 20, year == "1995"),
                aes(x =  bin_midpoint, y = abundance_ratio_pioneerslow, fill = abundance_ratio_pioneerslow),
                method = "lm", color = "black", alpha = 0.2 ) +
    geom_errorbar(data = abun_ratio_pioneerslow_range  %>% filter(bin_midpoint <= abun_max), 
                  aes(x = bin_midpoint, ymin = min, ymax = max), width = 0, color = "black") +
    geom_point(data = ratios %>% filter(min_n_individuals_pioneerslow >= 20, year == "1995"),
               aes(x =  bin_midpoint, y = abundance_ratio_pioneerslow, fill = abundance_ratio_pioneerslow),
               shape = 21, stroke = 0.5, size = 4, color = "black") +
    scale_slow_tall  +
    new_scale_fill() +
    stat_smooth(data = ratios %>% filter(min_n_individuals_fastslow >= 20, year == "1995"),
                aes(x =  bin_midpoint, y = abundance_ratio_fastslow, fill = abundance_ratio_fastslow),
                method = "lm", color = "black", alpha = 0.2 ) +
    geom_errorbar(data = abun_ratio_fastslow_range, aes(x = bin_midpoint, ymin = min, ymax = max), width = 0, color = "black") +
    geom_point(data = ratios %>% filter(min_n_individuals_fastslow >= 20, year == "1995"),
               aes(x =  bin_midpoint, y = abundance_ratio_fastslow, fill = abundance_ratio_fastslow),
               shape = 21, stroke = 0.5, size = 4, color = "black") +
    scale_slow_fast 
)


p1 <- set_panel_size(abun_ratio , width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


#----------------------------------------------------------------------------------
#--------------------------  Fig 2D:  Richness Ratio -----------------------------
#----------------------------------------------------------------------------------
ratio_pionslow_range <- ratios %>%
  filter(min_n_individuals_pioneerslow >= 20) %>%
  group_by(bin_midpoint) %>%
  summarize(
    min = min(richness_ratio_pioneerslow),
    max = max(richness_ratio_pioneerslow))

rich_ratio_fastslow_range <- ratios %>%
  filter(min_n_individuals_fastslow >= 20) %>%
  group_by(bin_midpoint) %>%
  summarize(
    min = min(richness_ratio_fastslow),
    max = max(richness_ratio_fastslow))

rich_data = ratios %>% filter(min_n_individuals_pioneerslow >= 20, min_n_individuals_fastslow >= 20, year == "1995")

rich_max = max(rich_data$bin_midpoint)

(rich_ratio  <- 
    ggplot() + 
    scale_x_log10(limits = c(.9, 100), breaks = c(1, 10, 100), 
                  position = "bottom", expression(paste('Stem Diameter (cm)'))) + 
    scale_y_log10(labels = signif,
                  breaks = c(0.01, 0.1, 1, 10), limits = c(0.03, 10),
                  position = "right", name = expression("Richness Ratio")) + 
    theme_plant() +
    stat_smooth(data = ratios %>% filter(min_n_individuals_pioneerslow >= 20, year == "1995"),
                aes(x =  bin_midpoint, y = richness_ratio_pioneerslow, fill = richness_ratio_pioneerslow),
                method = "lm", color = "black", alpha = 0.2 ) +
    geom_errorbar(data = ratio_pionslow_range %>% filter(bin_midpoint <= rich_max), aes(x = bin_midpoint, ymin = min, ymax = max), width = 0, color = "black") +
    geom_point(data = ratios %>% filter(min_n_individuals_pioneerslow >= 20, year == "1995"),
               aes(x =  bin_midpoint, y = richness_ratio_pioneerslow, fill = richness_ratio_pioneerslow),
               shape = 21, stroke = 0.5, size = 4, color = "black") +
    scale_slow_tall  +
    new_scale_fill() +
    stat_smooth(data = ratios %>% filter(min_n_individuals_fastslow >= 20, year == "1995"),
                aes(x =  bin_midpoint, y = richness_ratio_fastslow, fill = richness_ratio_fastslow),
                method = "lm", color = "black", alpha = 0.2 ) +
    geom_errorbar(data = rich_ratio_fastslow_range %>% filter(bin_midpoint <= rich_max), aes(x = bin_midpoint, ymin = min, ymax = max), width = 0, color = "black")  +
    geom_point(data = ratios %>% filter(min_n_individuals_fastslow >= 20, year == "1995"),
               aes(x =  bin_midpoint, y = richness_ratio_fastslow, fill = richness_ratio_fastslow),
               shape = 21, stroke = 0.5, size = 4, color = "black") +
    scale_slow_fast + theme_no_x()
)

p1 <- set_panel_size(rich_ratio, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


#----------------------------------------------------------------------------------
#--------------------------  Fig 2F: Productivity  Ratio --------------------------
#----------------------------------------------------------------------------------
#---- LL pioneer vs slow abundance ------
prod_ratio_pioneerslow_range <- ratios %>%
  group_by(bin_midpoint) %>%
  filter(min_n_individuals_pioneerslow >= 20) %>%
  summarize(
    min = min(production_ratio_pioneerslow),
    max = max(production_ratio_pioneerslow))

prod_ratio_fastslow_range <- ratios %>%
  group_by(bin_midpoint) %>%
  filter(min_n_individuals_fastslow >= 20) %>%
  summarize(
    min = min(production_ratio_fastslow),
    max = max(production_ratio_fastslow))

prod_data = ratios %>% filter(min_n_individuals_pioneerslow >= 20, min_n_individuals_fastslow >= 20, year == "1995")

prod_max = max(prod_data$bin_midpoint)

(prod_ratio  <-
    ggplot() + theme_plant() +
    scale_slow_tall  +
    scale_x_log10(limits = c(.9, 100), breaks = c(1, 10, 100), 
                  position = "bottom", 
                  expression(paste('Stem Diameter (cm)'))) + 
    scale_y_log10(labels = signif, 
                  breaks = c(0.01, 0.1, 1, 10), limits = c(0.03, 10),
                  position = "right", name = expression("Productivity Ratio")) + 
    stat_smooth(data = ratios %>%  filter(min_n_individuals_pioneerslow >= 20, year == "1995"),
                aes(x =  bin_midpoint, y = production_ratio_pioneerslow, fill = production_ratio_pioneerslow),
                method = "lm", color = "black", alpha = 0.2 ) +
    geom_errorbar(data = prod_ratio_pioneerslow_range %>% filter(bin_midpoint <= prod_max), aes(x = bin_midpoint , ymin = min, ymax = max), width = 0, color = "black")  +
    geom_point(data = ratios %>% filter(min_n_individuals_pioneerslow >= 20, year == "1995"),
               aes(x =  bin_midpoint, y = production_ratio_pioneerslow, fill = production_ratio_pioneerslow),
               shape = 21, stroke = 0.5, size = 4, color = "black") +
    new_scale_fill() +
    scale_slow_fast +
    stat_smooth(data = ratios %>% 
                  filter(min_n_individuals_fastslow >= 20, year == "1995"),
                aes(x =  bin_midpoint, y = production_ratio_fastslow, fill = production_ratio_fastslow),
                method = "lm", color = "black", alpha = 0.2 ) +
    geom_errorbar(data = prod_ratio_fastslow_range %>% filter(bin_midpoint <= prod_max), aes(x = bin_midpoint, ymin = min, ymax = max), width = 0, color = "black")  +
    geom_point(data = ratios %>% 
                 filter(min_n_individuals_fastslow >= 20, year == "1995"),
               aes(x =  bin_midpoint, y = production_ratio_fastslow, fill = production_ratio_fastslow),
               shape = 21, stroke = 0.5, size = 4, color = "black") 
)


p1 <- set_panel_size(prod_ratio, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


##########################################################################################
############################ Fig 3: Life History PCA Plots ######################################
##########################################################################################


fgbci$fg5 <- match(fgbci$fg5, c(2,3,1,4,5))

# Correct these
fgbci$PC_slow_to_fast <- -fgbci$X1new
fgbci$PC_breeder_to_pioneer <- fgbci$X2new

binomial <- paste(fgbci$genus, fgbci$species, sep = " ")
binomial
unique(binomial)
guild_lookup <- data.frame(fg = c('fg1','fg2','fg3','fg4','fg5','all','unclassified'), 
                           fg_name = c('Fast','Tall','Slow','Short','Medium','All','Unclassified'))

guild_fills_fg2 = c("1" = "#BFE046", "2" =  "#267038", "3" = "#27408b", "4" = "#87Cefa", "5" = "gray93")
geom_size = 3
geom_size = 4
unique(fgbci$fg5)


#----------------------------------------------------------------------------------
#--------------------------  Fig 3A: Life History PCA values ----------------------
#----------------------------------------------------------------------------------


p <- ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, fill = as.factor(fg5))) +
  geom_point(shape = 21, size = geom_size, stroke = 0.3, color = "black") + 
  labs(x = 'Survivorship–Growth Tradeoff', y = 'Stature—Recruitment Tradeoff') +
  theme_plant_small() + theme(aspect.ratio = 0.75) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,3)) +
  scale_x_continuous(limits = c(-6,7), breaks = seq(-6,6,3)) +
  scale_fill_manual(values = guild_fills_fg2) 
p

p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p2)


#----------------------------------------------------------------------------------
#--------------------------  Fig 3B: Mean PCA vs Stem Diameter --------------------
#----------------------------------------------------------------------------------

n_indiv = obs_totalprod %>% filter(fg == "all") %>%
  select(bin_midpoint, bin_count) %>%
  rename(n_individuals = bin_count)

fastslowscore_bin_bydiam <- fastslowscore_bin_bydiam_byyear %>%
  left_join(fastslow_stats_bydiam_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Fast-Slow")

breederscore_bin_bydiam <- breederscore_bin_bydiam_byyear %>%
  left_join(breeder_stats_bydiam_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = 'Short-Tall')

score_bin_bydiam <- rbind(fastslowscore_bin_bydiam, breederscore_bin_bydiam)

# Fast vs Slow Axis
PCA_diam_fs <- score_bin_bydiam  %>%
  filter(year == 1995, n_individuals >= 20, ID == "Fast-Slow") %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  theme_plant() +
  geom_errorbar() + 
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black",
             aes(fill = mean)) +
  scale_slow_fast2 +
  scale_y_continuous(limits = c(-1.5, 1.25),breaks = c(-1,0,1), name = 'PCA Score') +
  scale_x_log10(name = 'Stem Diameter (cm)', limits = c(.8, 100)) 
PCA_diam_fs

# LLP vs SLR axis
PCA_diam_st <-  score_bin_bydiam  %>%
  filter(year == 1995, n_individuals >= 20, ID == "Short-Tall") %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  theme_plant() +
  geom_errorbar() + 
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black",
             aes(fill = mean)) +
  scale_short_tall +
  scale_x_log10(name = 'Stem Diameter (cm)', limits = c(.8,100)) +
  scale_y_continuous(limits = c(-1.5, 1.25),breaks = c(-1,0,1), name = 'PCA Score') 
PCA_diam_st


#----------------------------------------------------------------------------------
#--------------------------  Fig 3C: Mean PCA vs Light Intensity ------------------
#----------------------------------------------------------------------------------
fastslowscore_bin_bylight <- fastslowscore_bin_bylight_byyear %>%
  left_join(fastslow_stats_bylight_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Fast-Slow")

breederscore_bin_bylight <- breederscore_bin_bylight_byyear  %>%
  left_join(breeder_stats_bylight_byyear %>% select(year, bin_midpoint, n_individuals)) %>%
  mutate(ID = "Short-Tall")

PCA_score_by_light <- rbind(fastslowscore_bin_bylight,breederscore_bin_bylight)

error_bar_width <- .15


# LLP vs SLR axis
p <- PCA_score_by_light %>%
  filter(year == 1995, n_individuals >= 20, ID == "Short-Tall") %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed")+ 
  geom_errorbar(width = error_bar_width) + theme_plant() +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black", aes(fill = mean)) +
  scale_short_tall +
  scale_x_log10(limits = c(1.8,450), breaks = c(3, 30, 300),
                name = expression(paste('Light Intensity (W m'^-2,')'))) + 
  scale_y_continuous(limits = c(-1.5, 1.25), breaks = c(-1, 0, 1), name = 'PCA Score') 
p

p1 <- set_panel_size(p, width=unit(11,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


# Fast vs Slow Axis
p <- PCA_score_by_light %>%
  filter(year == 1995, n_individuals >= 20, ID == "Fast-Slow") %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_errorbar(width = error_bar_width) + theme_plant() +
  geom_point(shape = 21, size = 4.5,  stroke = .5,  color = "black", aes(fill = mean)) +
  scale_slow_fast2 +
  scale_x_log10(limits = c(1.8,450), breaks = c(3, 30, 300),
                name = expression(paste('Light Intensity (W m'^-2,')'))) + 
  scale_y_continuous(limits = c(-1.5, 1.25), breaks = c(-1, 0, 1), name = 'PCA Score')
p

p1 <- set_panel_size(p, width=unit(11,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


#------------------------------------------------------------------------------------------------
#-------- Fig S2: Mean PCA with credible interval and range across census years ------------------
#------------------------------------------------------------------------------------------------

id_mid = c("ID", "bin_midpoint")

pca_range <- score_bin_bydiam  %>% filter( n_individuals >= 20) %>%
  group_by(across(id_mid)) %>%
  summarize(
    min = min(mean),
    max = max(mean)
  )

#--------------------------- Confidence intervals 
(pca_score_ci = 
    ggplot() + 
    theme_plant() +
    geom_line(data = score_bin_bydiam %>% filter(year == 1995, n_individuals >= 20), aes(x = bin_midpoint, y = mean,  color = ID)) +
    geom_point(data = score_bin_bydiam %>% filter(year == 1995, n_individuals >= 20), aes(x = bin_midpoint, y = mean,fill = ID),
               shape = 21, size = 4,  stroke = .5,  color = "black") +
    geom_errorbar(data = score_bin_bydiam %>% filter(year == 1995, n_individuals >= 20), aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max), width = 0, color = "forestgreen") +
    scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey")) +
    scale_color_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey50")) +
    scale_x_log10(limits = c(1, 100), breaks = c(1, 10, 100), 
                  name =NULL) + 
    scale_y_continuous(name = expression("PCA Score")) + theme_no_x() +
    theme_plant() )

p1 <- set_panel_size(pca_score_ci, width=unit(11,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

#---------------------- mean range over census years
(pca_score_range = 
    ggplot() + 
    theme_plant() +
    geom_line(data = score_bin_bydiam %>% filter(year == 1995, n_individuals >= 20), aes(x = bin_midpoint, y = mean,  color = ID)) +
    geom_point(data = score_bin_bydiam %>% filter(year == 1995, n_individuals >= 20), aes(x = bin_midpoint, y = mean,fill = ID),
               shape = 21, size = 4,  stroke = .5,  color = "black") +
    geom_errorbar(data = pca_range, aes(x = bin_midpoint, ymin = min, ymax = max), width = 0, color = "dodgerblue") +    
    scale_fill_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey")) +
    scale_color_manual(values = c("Short-Tall" = "black", "Fast-Slow" = "grey50")) +
    scale_x_log10(limits = c(1, 100), breaks = c(1, 10, 100), 
                  name = expression(paste('Diameter (cm)'))) + 
    scale_y_continuous(name = expression("PCA Score")) + 
    theme_plant() )


p1 <- set_panel_size(pca_score_range, width=unit(11,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)


######################################################################################################
#########################        Figure 4            #################################################
######################################################################################################


#------------------------------------------------------------------
#------------------ Fig 4A: Richness ~ Abundance ------------------
#------------------------------------------------------------------


max_size_fg <- obs_richnessbydiameter %>%
  filter(n_individuals > 20) %>%
  group_by(fg) %>%
  dplyr::summarize(max_size = max(abundance_by_bin_width)) 

min_size_fg <- obs_richnessbydiameter %>%
  filter(n_individuals > 20) %>%
  group_by(fg) %>%
  dplyr::summarize(min_size = min(abundance_by_bin_width)) 

fitted_richnessvsabundance_filt <- fitted_richnessvsabundance %>%
  left_join(max_size_fg, by = "fg") %>%
  left_join(min_size_fg, by = "fg") %>%
  group_by(fg) %>% 
  filter(abundance_by_bin_width >= min_size) %>%
  filter(abundance_by_bin_width <= max_size)

(rich_abun <- ggplot() + 
    theme_plant() +
    scale_x_log10(name = expression(paste("Abundance (cm"^-1,")")), 
                  breaks = c(1,  100, 10000), labels = c(1, 100, "10,000"), limit = c(1, 80000)) + 
    scale_y_log10(labels = signif, limit = c(0.3, 300),
                  position = "left", name = expression(paste("Richness (cm"^-1,")"))) + 
    scale_fill_manual(values = guild_fills_fg) +
    scale_color_manual(values = guild_colors_fg) +
    geom_jitter(data = obs_richnessbydiameter %>%   arrange(desc(fg)) %>% 
                  filter(n_individuals > 20) %>% filter(!fg %in% c('unclassified', 'all')),  
                aes(x = abundance_by_bin_width, y = richness_by_bin_width,
                    fill = fg, color = fg),
                shape = 21, size = 4, color = "black", width = 0)  +
    theme(axis.title.y = element_text(vjust = -3)) + 
    geom_ribbon(data = fitted_richnessvsabundance_filt %>% 
                  filter(!fg %in% c('unclassified', 'all')) %>% arrange(desc(fg)),
                aes(x = abundance_by_bin_width, ymin = q025, ymax = q975, fill = fg, color = NA), alpha = 0.2, col = NA) + 
    geom_line(data = fitted_richnessvsabundance_filt %>% 
                filter(!fg %in% c('unclassified', 'all')) %>%  arrange(desc(fg)),
              aes(x = abundance_by_bin_width, y = q50, color = fg), size = 0.5)
)


#------------------------------------------------------------------------------
#---------- Fig S5: Richness ~ Abundance in other census year ------------------
#-------------------------------------------------------------------------------

Year = 2000 #adjust for other years

(p <- 
    ggplot(data = binned_data %>% arrange(desc(fg)) %>% filter(abundance >= 20, year == Year, !fg %in% c('unclassified', 'all')), 
           aes(x = abundance_by_bin_width, y = richness_by_bin_width, fill = fg, color = fg)) + 
   theme_plant() +
   scale_x_log10(name = expression(paste("Abundance (cm"^-1,")")), 
                 breaks = c(1,  100, 10000), labels = c(1, 100, "10,000")) + 
   scale_y_log10(labels = signif,
                  name = expression(paste("Richness (cm"^-1,")"))) +
   scale_fill_manual(values = guild_fills_fg) +
   scale_color_manual(values = guild_colors_fg) +
   geom_point(shape = 21, size = 4, color = "black")  +
   theme(axis.title.y = element_text(vjust = -3)) + 
   stat_smooth(size = 0.5, alpha = 0.2, method = "lm"))


#----------------------------------------------------
#---------- Fig 4B Diameter Growth ------------------
#----------------------------------------------------

#used modified mass growth function
fg_mid = c("fg", "bin_midpoint")
diam_growth_range <- obs_indivdiamgrowth %>%
  group_by(across(fg_mid)) %>%
  filter(!fg %in% c("unclassified", "all"), mean_n_individuals >= 20) %>%
  summarize(
    min = min(mean),
    max = max(mean)
  )
diam_growth_range

fitted_diam_growth <- fitted_diamgrowthbydiameter

p <- plot_diam(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3','fg4','fg5'),
               model_fit = PROD,
               x_limits = c(.9, 160),
               x_breaks = c(1, 10, 100),
               y_limits = c(0.02, 1.3),
               y_breaks = c(0.03, 0.1, 0.3, 1),
               y_labels = c(0.03, 0.1, 0.3, 1),
               error_bar_width = 0,
               dodge_width = 0.05,
               obsdat = obs_indivdiamgrowth, 
               preddat = fitted_indivdiamgrowth,
               plot_abline = FALSE,
               x_name = 'Stem Diameter (cm)',
               y_name = expression(paste('Diameter Growth (cm yr'^-1,')')))
p

(p1 <- p + theme(axis.text.x = element_text(), axis.ticks.x = element_line()) + 
    annotation_custom(grob_fast) + annotation_custom(grob_tall) + annotation_custom(grob_medium) + 
    annotation_custom(grob_slow) + annotation_custom(grob_short) +
    labs(x = 'Stem Diameter (cm)') + theme(plot.margin = grid::unit(c(0,1,0,0), "cm")) +theme_plant_small())



#-----------------------------------------------------------------------------
#---------- Fig S4: Diameter Growth over other census years ------------------
#----------------------------------------------------------------------------

p <- plot_diam2(year_to_plot = 1995,
               fg_names = c('fg1','fg2','fg3','fg4','fg5'),
               model_fit = PROD,
               x_limits = c(.9, 160),
               x_breaks = c(1, 10, 100),
               y_limits = c(0.02, 1.3),
               y_breaks = c(0.03, 0.1, 0.3, 1),
               y_labels = c(0.03, 0.1, 0.3, 1),
               error_bar_width = 0,
               dodge_width = 0.05,
               obsdat = obs_indivdiamgrowth, 
               preddat = fitted_indivdiamgrowth,
               plot_abline = FALSE,
               x_name = 'Stem Diameter (cm)',
               y_name = expression(paste('Diameter Growth (cm yr'^-1,')')))
p

#---------------------------------------------------------------------------------------------------
#--------------------------   Fig 4C: Mortality   --------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# Read pre-created mortality bins for plotting to compare to fitted values
binned_data$mortality_by_bin_width <- binned_data$mortality/(binned_data$bin_max - binned_data$bin_min)
fitted_values <- fitted_mortalitybydiameter

# Filter the binned data to only have error bars for the years where we have a 1995 data point.
binned_data2 <- binned_data %>%
  filter(abundance > 20) %>%
  group_by(fg, bin) %>%
  mutate(has1995 = 1995 %in% year) %>%
  ungroup %>%
  filter(has1995) 


(p<- ggplot(binned_data2 %>% filter(fg %in% paste0('fg', 1:5), abundance > 20), aes(x=bin_midpoint, y=mortality, color = fg, fill = fg)) +
    theme_plant() + 
    scale_color_manual(values = guild_colors_fg) +
    scale_fill_manual(values = guild_fills_fg) +
    scale_y_continuous(trans = "logit", 
                       breaks = c(0.01, 0.03, 0.1,  0.3), 
                       limits = c(0.008, 0.3),
                       name = expression(paste("Mortality (yr"^-1,")"))) +
    scale_x_log10(name = parse(text = 'Stem~Diameter~(cm)'), 
                  breaks = c(1, 10, 100), 
                  limits = c(.9, 160)) +
    stat_summary(geom = 'errorbar', fun.max = max, fun.min = min, position = position_dodge(width= .05), width = 0) +
    geom_ribbon(data = fitted_values, aes(x = dbh, y = q50, ymin = q025, ymax = q975),  alpha = 0.4, color = NA) +
    geom_line(data = fitted_values, aes(x = dbh, y = q50, group = fg, color = fg)) +
    geom_point(aes(y = mortality, fill = fg), data = binned_data %>% filter(fg %in% paste0('fg', 1:5), abundance > 20, year == 1995), 
               position = position_dodge(width=.03), color = "black", shape = 21, size = 4) 
)

p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p2)  

#---------------------------------------------------------------------------------------------------
#--------------------------   Fig S4: Mortality over other census years  --------------------------
#---------------------------------------------------------------------------------------------------

(p<- ggplot(binned_data2 %>% filter(fg %in% paste0('fg', 1:5), abundance > 20), aes(x=bin_midpoint, y=mortality, color = fg, fill = fg)) +
   theme_plant() + 
   scale_color_manual(values = guild_colors_fg) +
   scale_fill_manual(values = guild_fills_fg) +
   scale_y_continuous(trans = "logit", 
                      breaks = c(0.01, 0.03, 0.1,  0.3), 
                      limits = c(0.008, 0.3),
                      name = expression(paste("Mortality (yr"^-1,")"))) +
   scale_x_log10(name = parse(text = 'Stem~Diameter~(cm)'), 
                 breaks = c(1, 10, 100), 
                 limits = c(.9, 160)) +
   stat_summary(geom = 'errorbar', fun.max = max, fun.min = min, position = position_dodge(width= .05), width=0) + 
   geom_line(data = fitted_values, aes(x = dbh, y = q50, group = fg, color = fg)) 
)

p2 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p2)  

#---------------------------------------------------------------------------------------------------
#--------------------------   Fig 4D: Mass Growth  --------------------------------------------------------
#---------------------------------------------------------------------------------------------------

mass_growth_range <- obs_indivprod %>%
  filter(!fg %in% c("unclassified", "all"), mean_n_individuals >= 20) %>%
  group_by(across(fg_mid )) %>%
  summarize(
    min = min(mean),
    max = max(mean)
  )
mass_growth_range


obs_indivprod <- obs_indivprod %>%
  filter(mean_n_individuals >= 20)

(p <- plot_growth(year_to_plot = 1995,
                  fg_names = c('fg1','fg2','fg3','fg4','fg5'),
                  model_fit = PROD,
                  x_limits = c(0.9, 185),
                  y_limits = c(0.003, 200),
                  y_breaks = c(0.01, 0.1, 1, 10, 100),
                  plot_errorbar = F,
                  error_min = 'q25',
                  error_max = 'q75',
                  error_bar_width = 0,
                  y_labels = c( 0.01, 0.1, 1, 10, 100),
                  abline_slope = 2, 
                  abline_intercept = 10, 
                  dodge_width = 0.05))


p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)

#---------------------------------------------------------------------------------------------------
#--------------------------   Fig S4: Mass Growth over other census years -------------------------
#---------------------------------------------------------------------------------------------------

(p <- plot_growth2(year_to_plot = 1995,
                  fg_names = c('fg1','fg2','fg3','fg4','fg5'),
                  model_fit = PROD,
                  x_limits = c(0.9, 185),
                  y_limits = c(0.003, 200),
                  y_breaks = c(0.01, 0.1, 1, 10, 100),
                  plot_errorbar = F,
                  error_min = 'q25',
                  error_max = 'q75',
                  error_bar_width = 0,
                  y_labels = c( 0.01, 0.1, 1, 10, 100),
                  abline_slope = 2, 
                  abline_intercept = 10, 
                  dodge_width = 0.05))


p1 <- set_panel_size(p, width=unit(10.25,"cm"), height=unit(8,"cm"))
grid.newpage()
grid.draw(p1)



##################################################################################################
########################################   Figure 5 ##############################################
##################################################################################################


#---------------------------------------------------------------------------------------------------
#--------------------------   Fig 5A: Intercepted light -------------------------
#---------------------------------------------------------------------------------------------------

grob_text <- grobTree(textGrob("Solar Equivalence", x = 0.27, y = 0.87, hjust = 0,
                               gp = gpar(col = "gold3", fontsize = 18))) 

grob_a <- grobTree(textGrob("a", x = 0.06, y = 0.91, gp = gpar(col = "black", fontsize = 25, fontface = "bold")))
totallightbins_fg <- totallightbins_fg %>%
  filter(bin_count >= 20)


tot_light <- plot_totalprod(year_to_plot = 1995,
                            fg_names = c('fg1','fg2','fg3', 'fg4', 'fg5', 'all'),
                            model_fit_density = DENS, 
                            model_fit_production = PROD,
                            x_limits = c(0.8, 200),
                            y_limits = c(100, 200000),
                            geom_size = 3.5,
                            y_breaks = c(100, 1000, 10000, 100000),
                            x_breaks = c(1, 10, 100),
                            y_labels = c("0.1", "1", "10", "100"),
                            y_name = expression(paste('Total Light Intercepted (kW cm'^-1,' ha'^-1,')')),
                            preddat = fitted_totallight,
                            obsdat = totallightbins_fg,
                            plot_abline = FALSE)


tot_light2 <- tot_light  + 
  scale_y_continuous(position = "right", trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("0.1", "1", "10", "100"), 
                     limits = c(200, 200000),
                     name = expression(atop('Total Light Intercepted',paste('(kW cm'^-1,' ha'^-1,')'))))  +
  scale_x_log10(name = "Stem Diameter (cm)", limits = c(0.8, 200), position = "top") +
  theme(aspect.ratio = 0.75) + 
  geom_abline(intercept = log10(70000), slope = 0, color = "#C9B074",
              linetype = "dashed", size = 0.75) +
  annotation_custom(grob_text)
plot(tot_light2)
p_tot_light <- tot_light2


#---------------------------------------------------------------------------------------------------
#--------------------------   Fig 5B: Comparisons of slopes -------------------------
#---------------------------------------------------------------------------------------------------
xvalues <- params %>% 
  filter(variable == 'density', model == DENS) %>%
  group_by(fg) %>%
  summarize(mid_cutoff = (mean[parameter == 'tau_high'] + mean[parameter == 'tau_low'])/2)


growth_slopes_atmiddle <- growth_slopes %>% 
  filter(variable == 'total_production', dens_model == DENS, prod_model == PROD) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

light_slopes_atmiddle <- light_slopes %>% 
  filter(variable == 'total_incoming_light', dens_model == DENS, prod_model == PROD) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

fg_full_names <- c('Fast', 'Tall', 'Slow', 'Short', 'Medium', 'All', 'Unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

allslopes <- rbind(growth_slopes_atmiddle, light_slopes_atmiddle) %>%
  ungroup %>%
  mutate(fg = factor(fg, levels = fgs, labels = fg_full_names))
grobb <- grobTree(textGrob("b", x = 0.04, y = 0.9,  hjust = 0,
                           gp = gpar(col = "black", fontsize = 25, fontface = "bold"))) 
grob1 <- grobTree(textGrob("Solar", x = 0.68, y = 0.94, hjust = 0,
                           gp = gpar(col = "gold3", fontsize = 18))) 
grob2 <- grobTree(textGrob("Production", x = 0.68, y = 0.86, hjust = 0,
                           gp = gpar(col = "darkgreen", fontsize = 18)))
grob3 <- grobTree(textGrob("Energy Equivalence", x = 0.24, y = 0.51, hjust = 0,
                           gp = gpar(col = "black", fontsize = 18))) 



# compare slopes
slopes <- ggplot(allslopes %>% filter(!fg %in% 'Unclassified'), 
                 aes(x = fg, y = q50, ymin = q025, ymax = q975, fill =  variable, color =variable)) +
  theme_plant() +
  geom_hline(yintercept = 0, linetype = 'dashed', size = .75) +
  geom_errorbar(position = position_dodge(width = 0.6), size = 0.75, width = 0) +
  geom_point(position = position_dodge(width = 0.6), shape = 21, size = 3.5, color = "black", stroke = 0.5) +
  geom_errorbar(data = allslopes %>% filter(fg %in% c('Tall', 'Slow', 'All')), 
                position = position_dodge(width = 0.6), size = 0.75, width = 0) +
  scale_y_continuous(position = "right", ) +
  labs(x = NULL) +
  scale_y_continuous(name = "Scaling Slope", position = "right", limits = c(-1.05, 1.3)) +
  scale_fill_manual(values = c('gold1', 'darkolivegreen3')) +
  scale_color_manual(values = c('gold3', 'darkgreen')) +
  theme(axis.text.x = element_text(angle = 25,  vjust = .7, face = "italic", size = 18)) +
  annotation_custom(grob1) + annotation_custom(grob2) + annotation_custom(grob3) 
slopes

#combine figures
g_slopes <- ggplotGrob(slopes)
g_tot_light <- ggplotGrob(p_tot_light)
combo <- rbind(g_tot_light, g_slopes, size = "first")
combo$widths <- unit.pmax(g_tot_light$widths,g_slopes$widths)
grid.newpage()
grid.draw(combo)


#---------------------------------------------------------------------------------------------------
#--------------------------   Fig S6: Light Capture Scaling -------------------------
#---------------------------------------------------------------------------------------------------

hex_scale_log_colors <- scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), bias = 1)(50),
                                             trans = 'log', name = 'Individuals', breaks = c(1, 10, 100, 1000, 10000), 
                                             labels = c(1, 10 , 100, 1000, 10000), limits=c(1, 10000))

exl <- expression(atop('Intercepted Light', paste('per Individual (W)')))
exd <- 'Diameter (cm)'

labels = trans_format("log10", math_format(10^.x))

indiv_light <- ggplot() +
  theme_plant() +
  scale_x_log10(name = exd, limits = c(.9, 200)) +
  scale_y_log10(name = exl, breaks = c(1,100,10000, 1000000), limits = c(1,1000000), 
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hex(data = alltree_light_95, aes(x = dbh_corr, y = light_received)) +
  geom_pointrange(data = unscaledlightbydbhcloudbin_fg %>% filter(fg %in% 'all'), 
                  aes(x = dbh_bin, y = mean, ymin = mean, ymax = mean), size = .5) +
  theme(legend.position = "right", legend.text = element_text(size = 15), legend.title = element_text(size = 16))+
  hex_scale_log_colors +
  geom_smooth(data = alltree_light_95, aes(x = dbh_corr, y = light_received),
              method = "lm", color = "black", size = .7) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))
indiv_light 


