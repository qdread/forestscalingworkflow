
#----------Scaling Density Plot Function ----------

plot_dens <- function(year_to_plot = 1995, 
                      fg_names = c("fg1", "fg2", "fg3", "fg4", "fg5", "all"),
                      model_fit = 1, 
                      x_limits, 
                      x_breaks = c(1, 3, 10, 30, 100, 300),
                      y_limits, 
                      y_breaks, 
                      y_labels, 
                      fill_names = guild_fills_all, #c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                      color_names = guild_colors_all, #c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),  
                      x_name = "Stem Diameter (cm)",
                      y_name = expression(paste("Density (n ha"^-1, "cm"^-1, ")")), 
                      geom_size = 4, 
                      obsdat = obs_dens, 
                      preddat = pred_dens, 
                      plot_abline = TRUE, 
                      abline_slope = -2, 
                      abline_intercept = -10,
                      dodge_width = 0.03) 
{
  pos <-  ggplot2::position_dodge(width = dodge_width)
  obsdat <- obsdat %>% dplyr::filter(fg %in% fg_names, year == year_to_plot, bin_count >= 20) %>% 
    dplyr::filter(bin_value > 0) %>% 
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1')))
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(bin_midpoint), 
                     max_obs = max(bin_midpoint))
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% 
    dplyr::filter(dens_model %in% model_fit, fg %in% fg_names, 
                  year == year_to_plot) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>%  
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1')))
  p <- ggplot2::ggplot() + 
    ggplot2::geom_ribbon(data = preddat, 
                         ggplot2::aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), 
                         alpha = 0.4)
  if (plot_abline) {
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", linetype = "dashed", 
                                  size = 0.75)
  }
  p + ggplot2::geom_line(data = preddat, 
                         ggplot2::aes(x = dbh, y = q50, group = fg, color = fg)) + 
    ggplot2::geom_line(data = preddat[preddat$fg == "fg5", ], 
                       ggplot2::aes(x = dbh, y = q50), color = "gray") + 
    ggplot2::geom_errorbar(data = dens_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), width = 0, position = pos) +
    ggplot2::geom_point(data = obsdat, #position = pos,
                        ggplot2::aes(x = bin_midpoint, y = bin_value, group = fg, fill = fg), 
                        size = geom_size, shape = 21, color = "black") + 
    # ggplot2::geom_errorbar(data = dens_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), width = 0, position = pos) +
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) + 
    ggplot2::scale_color_manual(values = color_names) + 
    ggplot2::scale_fill_manual(values = fill_names) + 
    theme_plant() + theme_no_x()
}

#-------------- Supplemental Density Function for Ranges ---------------------------
plot_dens2 <- function(year_to_plot = 1995, 
                       fg_names = c("fg1", "fg2", "fg3", "fg4", "fg5", "all"),
                       model_fit = 1, 
                       x_limits, 
                       x_breaks = c(1, 3, 10, 30, 100, 300),
                       y_limits, 
                       y_breaks, 
                       y_labels, 
                       fill_names = guild_fills_all, #c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                       color_names = guild_colors_all, #c("black","#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),  
                       x_name = "Stem Diameter (cm)",
                       y_name = expression(paste("Density (n ha"^-1, "cm"^-1, ")")), 
                       geom_size = 4, 
                       obsdat = obs_dens, 
                       preddat = pred_dens, 
                       plot_abline = TRUE, 
                       abline_slope = -2, 
                       abline_intercept = -10,
                       dodge_width = 0.05) 
{
  pos <-  ggplot2::position_dodge(width = dodge_width)
  obsdat <- obsdat %>% dplyr::filter(fg %in% fg_names, year == year_to_plot, bin_count >= 20) %>% 
    dplyr::filter(bin_value > 0) %>% 
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1')))
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(bin_midpoint), 
                     max_obs = max(bin_midpoint))
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% 
    dplyr::filter(dens_model %in% model_fit, fg %in% fg_names, 
                  year == year_to_plot) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>%  
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1')))
  p <- ggplot2::ggplot() + 
    ggplot2::geom_ribbon(data = preddat, 
                         ggplot2::aes(x = dbh, ymin = q025, ymax = q975, group = fg, fill = fg), 
                         alpha = 0)
  if (plot_abline) {
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", linetype = "dashed", 
                                  size = 0.75)
  }
  p + ggplot2::geom_line(data = preddat, 
                         ggplot2::aes(x = dbh, y = q50, group = fg, color = fg), size = 0.2) + 
    ggplot2::geom_line(data = preddat[preddat$fg == "fg5", ], 
                       ggplot2::aes(x = dbh, y = q50), color = "gray", size = 0.2) + 
    ggplot2::geom_errorbar(data = dens_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), size = .5, width = 0, position = pos) +
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) + 
    ggplot2::scale_color_manual(values = color_names) + 
    ggplot2::scale_fill_manual(values = fill_names) + 
    theme_plant() + theme_no_x()
}


#------------------ Productivity Function 
# Modifications:  changed plotting order, geom_size, geom colors
plot_totalprod <-function(year_to_plot = 1995, 
                           fg_names = c("fg1", "fg2", "fg3","fg4", "fg5", "all"), 
                           model_fit_density = 1, 
                           model_fit_production = 1, 
                           x_limits, 
                           x_breaks = c(1, 3, 10, 30, 100, 300), 
                           y_limits = c(0.03, 100), 
                           y_breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                           y_labels, 
                          fill_names = guild_fills_all, # c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                          color_names = guild_colors_all, #c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray"), 
                           x_name = "Diameter (cm)", 
                           y_name = expression(paste("Production (kg cm"^-1, " ha"^-1, "  yr"^-1, ")")),
                           geom_size = 4, 
                           obsdat = obs_totalprod, 
                           preddat = fitted_totalprod, 
                           plot_abline = TRUE, 
                           abline_slope = 0, 
                           dodge_width = 0.0,
                           abline_intercept = 2) 
{
  pos <-  ggplot2::position_dodge(width = dodge_width)
  obsdat <- obsdat %>% 
    dplyr::filter(fg %in% fg_names, year == year_to_plot, bin_count >= 20) %>% 
    dplyr::filter(bin_value > 0) %>% 
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1')))
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(bin_midpoint), 
                     max_obs = max(bin_midpoint)) 
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% 
    dplyr::filter(dens_model %in% model_fit_density, 
                  prod_model %in% model_fit_production, 
                  fg %in% fg_names, year == year_to_plot) %>% 
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs)
  
  p <- ggplot2::ggplot() + 
    ggplot2::geom_ribbon(data = preddat, 
                         aes(x = dbh, ymin = q025, ymax = q975, 
                             group = fg, fill = fg), alpha = 0.4, show.legend = F) + 
    ggplot2::geom_line(data = preddat, show.legend = F,
                       aes(x = dbh, y = q50, group = fg, color = fg)) + 
    ggplot2::geom_point(data = obsdat, position = pos,
                        aes(x = bin_midpoint, y = bin_value, group = fg, fill = fg), 
                        size = geom_size, color = "black", shape = 21) + 
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, 
                           limits = y_limits, breaks = y_breaks, labels = y_labels,  position = "right") + 
    ggplot2::scale_color_manual(values = guild_colors_all) + 
    ggplot2::scale_fill_manual(values = guild_fills_all) + 
    theme_plant() + theme_no_x() +
    theme(aspect.ratio = 1)
  if (plot_abline) 
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", 
                                  linetype = "dashed", size = 0.75)
  p
} 


#--------------supplemental just ranges ---------
plot_totalprod2 <-function(year_to_plot = 1995, 
                           fg_names = c("fg1", "fg2", "fg3","fg4", "fg5", "all"), 
                           model_fit_density = 1, 
                           model_fit_production = 1, 
                           x_limits, 
                           x_breaks = c(1, 3, 10, 30, 100, 300), 
                           y_limits = c(0.03, 100), 
                           y_breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                           y_labels, 
                           fill_names = guild_fills_all, # c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                           color_names = guild_colors_all, #c("black", "#BFE046", "#267038", "#27408b", "#87Cefa", "gray"), 
                           x_name = "Stem Diameter (cm)", 
                           y_name = expression(paste("Productivity (kg cm"^-1, " ha"^-1, "  yr"^-1, ")")),
                           geom_size = 4, 
                           obsdat = obs_totalprod, 
                           preddat = fitted_totalprod, 
                           plot_abline = TRUE, 
                           abline_slope = 0, 
                           dodge_width = 0.03,
                           abline_intercept = 20) 
{
  pos <-  ggplot2::position_dodge(width = dodge_width)
  obsdat <- obsdat %>% 
    dplyr::filter(fg %in% fg_names, year == year_to_plot, bin_count >= 20) %>% 
    dplyr::filter(bin_value > 0) %>% 
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1')))
  obs_limits <- obsdat %>% dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(bin_midpoint), 
                     max_obs = max(bin_midpoint)) 
  preddat <- preddat %>% dplyr::left_join(obs_limits) %>% 
    dplyr::filter(dens_model %in% model_fit_density, 
                  prod_model %in% model_fit_production, 
                  fg %in% fg_names, year == year_to_plot) %>% 
    arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs)
  
  p <- ggplot2::ggplot() + 
    ggplot2::geom_ribbon(data = preddat, 
                         aes(x = dbh, ymin = q025, ymax = q975, 
                             group = fg, fill = fg), alpha = 0, show.legend = F) + 
    ggplot2::geom_line(data = preddat, show.legend = F,
                       aes(x = dbh, y = q50, group = fg, color = fg), size = 0.5) + 
    ggplot2::geom_errorbar(position = position_dodge(width = dodge_width), data = prod_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg), width = 0) +
    #ggplot2::geom_point(data = obsdat, position = pos,
    #                   aes(x = bin_midpoint, y = bin_value, group = fg, fill = fg), 
    #                  size = geom_size, color = "black", shape = 21) + 
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, 
                           limits = y_limits, breaks = y_breaks, labels = y_labels,  position = "left") + 
    ggplot2::scale_color_manual(values = guild_colors_all) + 
    ggplot2::scale_fill_manual(values = guild_fills_all) + 
    theme_plant() + #theme_no_x() +
    theme(aspect.ratio = 1)
  if (plot_abline) 
    p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                  slope = abline_slope, color = "gray72", 
                                  linetype = "dashed", size = 0.75)
  p
} 

#---------------- plot diameter growth ---------------
plot_diam <- function (year_to_plot = 1995, 
                        fg_names = c("Fast", "Tall", "Slow",  "Short", "Medium"), 
                        model_fit = 1, 
                        x_limits, 
                        x_breaks = c(1, 3, 10, 30, 100, 300), 
                        y_limits, 
                        y_labels, 
                        y_breaks, 
                        fill_names = guild_fills_fg, # c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                        color_names = guild_colors_fg, #c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),
                        x_name = "Diameter (cm)", 
                        y_name = expression(paste("Diameter Growth (cm yr"^-1,")")),
                        average = "mean", 
                        plot_errorbar = FALSE, 
                        error_min = "ci_min",
                        error_max = "ci_max", 
                        error_bar_width = 0.1,
                        error_bar_thickness = 0.5, 
                        dodge_width = 0.03, 
                        dodge_errorbar = TRUE, 
                        geom_size = 4, 
                        obsdat = obs_indivprod, 
                        preddat = fitted_diam_growth, 
                        plot_abline = TRUE, 
                        abline_slope = 2, 
                        abline_intercept = -1.25) {
  pos <- if (dodge_errorbar) 
    ggplot2::position_dodge(width = dodge_width)
  else "identity"
  obsdat <- obsdat %>% 
    dplyr::filter(fg %in% fg_names, year ==  year_to_plot, 
                  !is.na(mean), mean_n_individuals >= 20) %>% 
    dplyr::group_by(bin_midpoint) %>% 
    dplyr::mutate(width = error_bar_width *  dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    arrange(desc(fg))
  obs_limits <- obsdat %>% 
    dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(bin_midpoint), 
                     max_obs = max(bin_midpoint))
  preddat <- preddat %>% 
    dplyr::left_join(obs_limits) %>% 
    dplyr::filter(prod_model %in% model_fit, 
                  fg %in% fg_names, year == year_to_plot) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>% 
    arrange(desc(fg))
  p <- ggplot2::ggplot() + 
     ggplot2::geom_ribbon(data = fitted_diam_growth  %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))), 
                         ggplot2::aes(x = dbh, ymin = q025, ymax = q975, 
                                     group = fg, fill = fg), alpha = 0.4) + 
    ggplot2::geom_line(data = fitted_diam_growth  %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))), 
                       ggplot2::aes(x = dbh, y = q50, group = fg, color = fg))
  if (plot_errorbar) {
    p <- p + ggplot2::geom_errorbar(data = obsdat %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))), 
                                    ggplot2::aes_string(x = "bin_midpoint", 
                                                        ymin = error_min, 
                                                        ymax = error_max, 
                                                        group = "fg", 
                                                        color = "fg", 
                                                        width = "width"), 
                                    position = pos, 
                                    size = error_bar_thickness)
  }
  p <- p + ggplot2::geom_line(data = fitted_diam_growth[fitted_diam_growth$fg == "Medium", ], 
                              ggplot2::aes(x = dbh, y = q50), color = "gray") + 
    ggplot2::geom_errorbar(data = diam_growth_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg), width = 0, position = pos) +
     ggplot2::geom_point(data = obsdat %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))),
                        ggplot2::aes_string(x = "bin_midpoint", 
                                            y = average, group = "fg", fill = "fg"), 
                       size = geom_size, color = "black", shape = 21, position = pos) + 
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) + 
    theme_no_x() + 
    ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")) + 
    ggplot2::scale_color_manual(values = color_names) + 
    ggplot2::scale_fill_manual(values = fill_names) + 
    theme_plant() +
    if (plot_abline) {
      p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                    slope = abline_slope, color = "gray72", linetype = "dashed", 
                                    size = 0.75)
    }
  p
}

#-----plot diameter growth ranges----
plot_diam2 <- function (year_to_plot = 1995, 
                       fg_names = c("Fast", "Tall", "Slow",  "Short", "Medium"), 
                       model_fit = 1, 
                       x_limits, 
                       x_breaks = c(1, 3, 10, 30, 100, 300), 
                       y_limits, 
                       y_labels, 
                       y_breaks, 
                       fill_names = guild_fills_fg, # c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
                       color_names = guild_colors_fg, #c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),
                       x_name = "Diameter (cm)", 
                       y_name = expression(paste("Diameter Growth (cm yr"^-1,")")),
                       average = "mean", 
                       plot_errorbar = FALSE, 
                       error_min = "ci_min",
                       error_max = "ci_max", 
                       error_bar_width = 0.1,
                       error_bar_thickness = 0.5, 
                       dodge_width = 0.03, 
                       dodge_errorbar = TRUE, 
                       geom_size = 4, 
                       obsdat = obs_indivprod, 
                       preddat = fitted_diam_growth, 
                       plot_abline = TRUE, 
                       abline_slope = 2, 
                       abline_intercept = -1.25) {
  pos <- if (dodge_errorbar) 
    ggplot2::position_dodge(width = dodge_width)
  else "identity"
  obsdat <- obsdat %>% 
    dplyr::filter(fg %in% fg_names, year ==  year_to_plot, 
                  !is.na(mean), mean_n_individuals >= 20) %>% 
    dplyr::group_by(bin_midpoint) %>% 
    dplyr::mutate(width = error_bar_width *  dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    arrange(desc(fg))
  obs_limits <- obsdat %>% 
    dplyr::group_by(fg) %>% 
    dplyr::summarize(min_obs = min(bin_midpoint), 
                     max_obs = max(bin_midpoint))
  preddat <- preddat %>% 
    dplyr::left_join(obs_limits) %>% 
    dplyr::filter(prod_model %in% model_fit, 
                  fg %in% fg_names, year == year_to_plot) %>% 
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                     dplyr::all_vars(. > min(y_limits))) %>% 
    dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>% 
    arrange(desc(fg))
  p <- ggplot2::ggplot() + 
    # ggplot2::geom_ribbon(data = fitted_diam_growth  %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))), 
    #                     ggplot2::aes(x = dbh, ymin = q025, ymax = q975, 
    #                                 group = fg, fill = fg), alpha = 0.4) + 
    ggplot2::geom_line(data = fitted_diam_growth  %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))), 
                       ggplot2::aes(x = dbh, y = q50, group = fg, color = fg))
  if (plot_errorbar) {
    p <- p + ggplot2::geom_errorbar(data = obsdat %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))), 
                                    ggplot2::aes_string(x = "bin_midpoint", 
                                                        ymin = error_min, 
                                                        ymax = error_max, 
                                                        group = "fg", 
                                                        color = "fg", 
                                                        width = "width"), 
                                    position = pos, 
                                    size = error_bar_thickness)
  }
  p <- p + ggplot2::geom_line(data = fitted_diam_growth[fitted_diam_growth$fg == "Medium", ], 
                              ggplot2::aes(x = dbh, y = q50), color = "gray") + 
    ggplot2::geom_errorbar(data = diam_growth_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg), width = 0, position = pos) +
    # ggplot2::geom_point(data = obsdat %>% arrange(factor(fg, levels = c('fg5','fg4','fg3','fg2','fg1'))),
    #                    ggplot2::aes_string(x = "bin_midpoint", 
    #                                        y = average, group = "fg", fill = "fg"), 
    #                   size = geom_size, color = "black", shape = 21, position = pos) + 
    ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
    ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels) + 
    theme_no_x() + 
    ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")) + 
    ggplot2::scale_color_manual(values = color_names) + 
    ggplot2::scale_fill_manual(values = fill_names) + 
    theme_plant() +
    if (plot_abline) {
      p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                    slope = abline_slope, color = "gray72", linetype = "dashed", 
                                    size = 0.75)
    }
  p
}
#------------ Plot individual mass growth -----------------
plot_growth <- 
  function (year_to_plot = 1995, 
            fg_names = c("Fast", "Tall", "Slow",  "Short", "Medium", "All"), 
            model_fit = 1, 
            x_limits, 
            x_breaks = c(1, 10, 100), 
            y_limits, 
            y_labels, 
            y_breaks, 
            fill_names = guild_fills_fg, # c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
            color_names = guild_colors_fg, #c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),
            x_name = "Stem Diameter (cm)", 
            y_name = expression(paste("Mass Growth (kg yr"^-1, ")")),
            average = "mean", 
            position = "right",
            plot_errorbar = FALSE, 
            error_min = "ci_min",
            error_max = "ci_max", 
            error_bar_width = 0.1,
            error_bar_thickness = 0.5, 
            dodge_width = 0.03, 
            dodge_errorbar = TRUE, 
            geom_size = 4, 
            obsdat = obs_indivprod, 
            preddat = fitted_indivprod, 
            plot_abline = TRUE, 
            abline_slope = 2, 
            abline_intercept = -1.25) {
    pos <- if (dodge_errorbar) 
      ggplot2::position_dodge(width = dodge_width)
    else "identity"
    obsdat <- obsdat %>% 
      dplyr::filter(fg %in% fg_names, year ==  year_to_plot, 
                    !is.na(mean), mean_n_individuals >= 20) %>% 
      dplyr::group_by(bin_midpoint) %>% 
      dplyr::mutate(width = error_bar_width *  dplyr::n()) %>% 
      dplyr::ungroup() %>% 
      arrange(desc(fg))
    obs_limits <- obsdat %>% 
      dplyr::group_by(fg) %>% 
      dplyr::summarize(min_obs = min(bin_midpoint), 
                       max_obs = max(bin_midpoint))
    preddat <- preddat %>% 
      dplyr::left_join(obs_limits) %>% 
      dplyr::filter(prod_model %in% model_fit, 
                    fg %in% fg_names, year == year_to_plot) %>% 
      dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                       dplyr::all_vars(. > min(y_limits))) %>% 
      dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>% 
      arrange(desc(fg))
    p <- ggplot2::ggplot() + 
      ggplot2::geom_ribbon(data = preddat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), #add hashtag to get supplemental plot 
                           ggplot2::aes(x = dbh, ymin = q025, ymax = q975, 
                                        group = fg, fill = fg), alpha = 0.4) + 
      ggplot2::geom_line(data = preddat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), 
                         ggplot2::aes(x = dbh, y = q50, group = fg, color = fg))
    if (plot_errorbar) {
      p <- p + ggplot2::geom_errorbar(data = obsdat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), 
                                      ggplot2::aes_string(x = "bin_midpoint", 
                                                          ymin = error_min, 
                                                          ymax = error_max, 
                                                          group = "fg", 
                                                          color = "fg", 
                                                          width = "width"), 
                                      position = pos, 
                                      size = error_bar_thickness)
    }
    p <- p + ggplot2::geom_line(data = preddat[preddat$fg == "Medium", ], 
                                ggplot2::aes(x = dbh, y = q50), color = "gray") + 
      ggplot2::geom_errorbar(data = mass_growth_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), position = pos, width = 0) +
      ggplot2::geom_point(data = obsdat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), #add hashtag to get supplemental plot 
                          ggplot2::aes_string(x = "bin_midpoint", 
                                              y = average, group = "fg", fill = "fg"), 
                          size = geom_size, color = "black", shape = 21) + # position = pos
      ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
      ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels, position = position) + 
      #theme_no_x() + 
      ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")) + 
      ggplot2::scale_color_manual(values = color_names) + 
      ggplot2::scale_fill_manual(values = fill_names) + 
      #theme_no_x() +
      theme_plant() #+ theme_no_x()
    if (plot_abline) {
      p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                    slope = abline_slope, color = "gray72", linetype = "dashed", 
                                    size = 0.75)
    }
    p
  }


#-------- plot growth ranges only -----
#------------ Plot individual mass growth -----------------
plot_growth2 <- 
  function (year_to_plot = 1995, 
            fg_names = c("Fast", "Tall", "Slow",  "Short", "Medium", "All"), 
            model_fit = 1, 
            x_limits, 
            x_breaks = c(1, 10, 100), 
            y_limits, 
            y_labels, 
            y_breaks, 
            fill_names = guild_fills_fg, # c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray87"), 
            color_names = guild_colors_fg, #c("#BFE046", "#267038", "#27408b", "#87Cefa", "gray"),
            x_name = "Stem Diameter (cm)", 
            y_name = expression(paste("Mass Growth (kg yr"^-1, ")")),
            average = "mean", 
            position = "right",
            plot_errorbar = FALSE, 
            error_min = "ci_min",
            error_max = "ci_max", 
            error_bar_width = 0.1,
            error_bar_thickness = 0.5, 
            dodge_width = 0.03, 
            dodge_errorbar = TRUE, 
            geom_size = 4, 
            obsdat = obs_indivprod, 
            preddat = fitted_indivprod, 
            plot_abline = TRUE, 
            abline_slope = 2, 
            abline_intercept = -1.25) {
    pos <- if (dodge_errorbar) 
      ggplot2::position_dodge(width = dodge_width)
    else "identity"
    obsdat <- obsdat %>% 
      dplyr::filter(fg %in% fg_names, year ==  year_to_plot, 
                    !is.na(mean), mean_n_individuals >= 20) %>% 
      dplyr::group_by(bin_midpoint) %>% 
      dplyr::mutate(width = error_bar_width *  dplyr::n()) %>% 
      dplyr::ungroup() %>% 
      arrange(desc(fg))
    obs_limits <- obsdat %>% 
      dplyr::group_by(fg) %>% 
      dplyr::summarize(min_obs = min(bin_midpoint), 
                       max_obs = max(bin_midpoint))
    preddat <- preddat %>% 
      dplyr::left_join(obs_limits) %>% 
      dplyr::filter(prod_model %in% model_fit, 
                    fg %in% fg_names, year == year_to_plot) %>% 
      dplyr::filter_at(dplyr::vars(dplyr::starts_with("q")), 
                       dplyr::all_vars(. > min(y_limits))) %>% 
      dplyr::filter(dbh >=  min_obs & dbh <= max_obs) %>% 
      arrange(desc(fg))
    p <- ggplot2::ggplot() + 
      #ggplot2::geom_ribbon(data = preddat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), #add hashtag to get supplemental plot 
       #                    ggplot2::aes(x = dbh, ymin = q025, ymax = q975, 
        #                                group = fg, fill = fg), alpha = 0.4) + 
      ggplot2::geom_line(data = preddat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), 
                         ggplot2::aes(x = dbh, y = q50, group = fg, color = fg))
    if (plot_errorbar) {
      p <- p + ggplot2::geom_errorbar(data = obsdat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), 
                                      ggplot2::aes_string(x = "bin_midpoint", 
                                                          ymin = error_min, 
                                                          ymax = error_max, 
                                                          group = "fg", 
                                                          color = "fg", 
                                                          width = "width"), 
                                      position = pos, 
                                      size = error_bar_thickness)
    }
    p <- p + ggplot2::geom_line(data = preddat[preddat$fg == "Medium", ], 
                                ggplot2::aes(x = dbh, y = q50), color = "gray") + 
      ggplot2::geom_errorbar(data = mass_growth_range, aes(x = bin_midpoint, ymin = min, ymax = max, color = fg ), position = pos, width = 0) +
      #ggplot2::geom_point(data = obsdat %>% arrange(factor(fg, levels = c('all', 'fg5','fg4','fg3','fg2','fg1'))), #add hashtag to get supplemental plot 
       #                   ggplot2::aes_string(x = "bin_midpoint", 
        #                                      y = average, group = "fg", fill = "fg"), 
         #                 size = geom_size, color = "black", shape = 21) + # position = pos
      ggplot2::scale_x_log10(name = x_name, limits = x_limits, breaks = x_breaks) + 
      ggplot2::scale_y_log10(name = y_name, limits = y_limits, breaks = y_breaks, labels = y_labels, position = position) + 
      #theme_no_x() + 
      ggplot2::theme(rect = ggplot2::element_rect(fill = "transparent")) + 
      ggplot2::scale_color_manual(values = color_names) + 
      ggplot2::scale_fill_manual(values = fill_names) + 
      #theme_no_x() +
      theme_plant() #+ theme_no_x()
    if (plot_abline) {
      p <- p + ggplot2::geom_abline(intercept = abline_intercept, 
                                    slope = abline_slope, color = "gray72", linetype = "dashed", 
                                    size = 0.75)
    }
    p
  }
