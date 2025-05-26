# Load libraries
library(tidyverse)
library(GenBinomApps) # to get Clopper-Pearson confidence interval function
library(cols4all)
library(future)
library(doFuture)
library(foreach)
library(progressr) # progress bars for long parallel computations
library(cowplot)
library(latex2exp)
library(fitdistrplus)
library(data.table)
library(multidplyr)
library(scales)
library(ggpubr)
library(egg)
library(ggh4x)
library(ggpointdensity)
library(FNN)

# Labeling functions
appender_R0 <- function(string) TeX(paste("$R_0 = $", string))
appender_sigma <- function(string) TeX(paste("$\\sigma = $", string))

# Reduce resolution of a vector by sub-sampling to length new_length
res_reduce <- function(df, new_length) {
  old_length <- length(unique(df))
  if (old_length < new_length) {
    warning("new_length is larger than old length. Vector will be unchanged")
    new_length <- old_length
  }
  ret <- unique(df)[seq(1, length(unique(df)), length.out = new_length)]
}

# Parameters ---- 
b <- 0.3 # biting rate
Tvh <- 0.5 # vector transmission-probability susceptible vector being infected after contact with infectious host
Nh <- 10000 # number of hosts in system
Nv <- 100000 # number of vectors in system
muv <- 0.1 # vector mortality rate
gammah <- 0.1 # host recovery rate

R0_from_Thv_function <- function(Thv) {
  sqrt(Thv * b^2 * Tvh * Nv / (Nh * gammah * muv))
}

# Sigma values (environmental noise level in 0-0.3) to test
sigmas <- seq(0, 1, by = 0.05)

# Final time point
max_time = 3650

# Load data ----
# all_df_modified = read_rds("./data/all_modified.rds")
# enviro_df_modified = read_rds("./data/enviro_modified.rds")

all_stats_df = read_rds("./data/all_stats.rds")
enviro_stats_df = read_rds("./data/enviro_stats.rds")
full_stats_df <- rbind(
  mutate(all_stats_df, type = "all"),
  mutate(enviro_stats_df, type = "enviro")
)
comp_stats_df = read_rds("./data/comp_stats.rds")
All_sims_plot_df = read_rds("./data/all_sims.rds")
comparison_trajectories = read_rds("./data/comp_trajectories.rds")

# Plots ----

R0s = c(0.95, 1.125, 2, 4.625)

R0_colors = c(
  # Hotter colors for R0 > 1
  rev(c4a("brewer.yl_or_br", sum(R0s > 1))),
  # Black for R0 == 1
  "black",
  # Cooler colors for R0 < 1
  c4a("scico.oslo", sum(R0s < 1))
)

nice_output_labeller = function(output_name) {
  output_label = case_when(
    output_name == "small_outbreak" ~ "Pr(Outbreak > 10 cases)",
    output_name == "big_outbreak" ~ "Pr(Outbreak > 100 cases)",
    output_name == "endemic" ~ "Pr(Transmission lasts 10 years)",
    output_name == "max_cases" ~ "Intensity",
    output_name == "duration" ~ "Outbreak duration",
    output_name == "max_time" ~ "Peak timing",
    output_name == "duration_dieout" ~ "Outbreak duration (non-endemic case)"
  )
  
  
  return(output_label)
}


# Generic plot function
generic_plot_function <- function(output_name, in_df) {
  
  R0_colors = rev(c4a("kovesi.bu_bk_rd", length(unique(in_df$R0_factor))))
  
  in_df %>% 
    ungroup() %>% 
    filter(name == output_name) %>% 
    # filter(R0_factor %in% c(0.95, 1.125, 2, 4.625)) %>% 
    ggplot(aes(x = sigma, y = mean, group = R0_factor)) +
    # Mean
    geom_line(aes(color = R0_factor), lwd = 1) +
    # Confidence intervals
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = R0_factor), alpha = 0.05) +
    scale_color_manual("R0",
                       # breaks = rev(R0s),
                       values = R0_colors) +
    scale_fill_manual("R0",
                      # breaks = rev(R0s),
                      values = R0_colors) +
    scale_x_continuous(TeX("Environmental noise strength [$\\sigma$]"),
                       limits = c(0,NA),
                       expand = c(0,0)) +
    scale_y_continuous(output_name,
                       limits = c(0,NA),
                       expand = c(0,0)) +
    theme_half_open()
}


# Create a "stretched" region at sigma = 0 to show what occurs with no environmental noise
# Create a "stretched" region at sigma = 0 to show what occurs with no environmental noise
stretch_sigma <- function(in_df, include_det = F) {
  # Step size
  sigma_step = (diff(unique(all_stats_df$sigma)))[1]
  
  # Filter values where sigma = 0
  zero_sigma <- in_df %>% filter(sigma == 0)
  
  # New sigma values
  negative_sigmas = -seq(0.0, 0.55, by = sigma_step)
  negative_sigmas = negative_sigmas[-length(negative_sigmas)]
  
  # Repeat the zero_sigma rows for each negative sigma value
  stretched <- zero_sigma[rep(1:nrow(zero_sigma), each = length(negative_sigmas)), ]
  
  # Assign the new sigma values
  stretched$sigma <- rep(negative_sigmas, times = nrow(zero_sigma))
  
  # Combine original data with stretched region
  out <- bind_rows(in_df, stretched)
  
  if (include_det) {
    det_zeros = in_df %>% filter(sigma == 0, type == "enviro")
    # New sigma values
    det_sigmas = -seq(0.5, 1.0, by = sigma_step)
    # det_sigmas = det_sigmas[-c(1:2)]
    
    # Repeat the zero_sigma rows for each negative sigma value
    det_stretched <- det_zeros[rep(1:nrow(det_zeros), each = length(det_sigmas)), ]
    
    # Assign the new sigma values
    det_stretched$sigma <- rep(det_sigmas, times = nrow(det_zeros))
    det_stretched = mutate(det_stretched, type = "all")
    
    out <- bind_rows(out, det_stretched)
  }
  return(out)
}


# Find a smoothed version of the absolute difference = 0 contour
smooth_zero_contour <- function(in_df, output_name, type_name, level_val) {
  
  if (type_name %in% c("all", "enviro")) {
    out_df = in_df %>% filter(type == type_name) %>% dplyr::select(-type)
  } else {
    out_df = in_df
  }
  var_name = colnames(out_df)[colnames(out_df) %in% c("mean", "abs_diff")]
  
  # Use nearest neighbor averaging to replace values of the differences
  
  smooth_level = 100
  nn_df = get.knnx(data = cbind(out_df$sigma, out_df$R0),
                   query = cbind(out_df$sigma, out_df$R0),
                   k = smooth_level)
  
  out_df$var_smooth <- rowMeans(matrix(out_df[[var_name]][nn_df$nn.index], ncol = smooth_level))
  
  sigma_vals = sort(unique(out_df$sigma))
  R0_vals = sort(unique(out_df$R0))
  var_mat = matrix(out_df$var_smooth, ncol = length(sigma_vals), byrow = TRUE)
  
  # Find contour for 0 from the averaged values
  contours <- contourLines(sigma_vals, R0_vals, var_mat, levels = level_val)
  
  # Convert contour output to dataframe for ggplot
  contour_df <- do.call(rbind, lapply(seq_along(contours), function(i) {
    data.frame(sigma = contours[[i]]$x, R0 = contours[[i]]$y, group = i)
  })) 
  
  if (!is.null(contour_df)) {
    
    if (type_name == "comp") {
      temp_df <- contour_df %>% filter(sigma == 0) %>% stretch_sigma(., include_det = F)
      temp_df <- rbind(temp_df, distinct(temp_df, R0, group) %>% mutate(sigma = -.5))
    }
    
    if (type_name == "all") {
      temp_df <- contour_df %>% filter(sigma == 0) %>% stretch_sigma(., include_det = F)
      temp_df <- rbind(temp_df, distinct(temp_df, R0, group) %>% mutate(sigma = -.5))
      temp_df_det <- in_df  %>% 
        filter(type == "enviro", sigma == 0, mean <= level_val)  %>% 
        stretch_sigma(., T) %>% 
        filter(sigma < -0.45, R0 == max(R0)) %>% 
        dplyr::select(sigma, R0) %>% 
        mutate(group = 1)
      temp_df_det <- rbind(temp_df_det, distinct(temp_df_det, R0, group) %>% mutate(sigma = -1.05))
      temp_df = rbind(temp_df_det, temp_df) %>% 
        arrange(sigma) 
    } 
    if (type_name == "enviro") {
      temp_df <- in_df  %>% 
        filter(type == "enviro", sigma == 0, mean <= level_val)  %>% 
        stretch_sigma(., T) %>% 
        filter(sigma > -0.51, R0 == max(R0)) %>% 
        dplyr::select(sigma, R0) %>% 
        mutate(group = 1) %>% arrange(sigma)
    }
    
    contour_df <- rbind(temp_df, contour_df)
  }
}

# Generic plot function
generic_heat_function <- function(output_name, type_name) {
  
  # Set heatmap palette to match unit of output
  palette = case_when(
    output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ "viridis",
    output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ "plasma",
    output_name %in% c("max_cases") ~ "turbo",
  )
  
  # Label legend appropriately
  leg_label = case_when(
    output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ "Probability",
    output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ "Time [years]",
    output_name %in% c("max_cases") ~ "Number of cases",
  )
  
  # Define zero proxy values
  zero_proxy = case_when(
    output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ 1E-5,
    output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ 1,
    output_name %in% c("max_cases") ~ 1,
  )
  
  noise_labels = ifelse(
    type_name == "all",
    list(c("Deterministic", "Demographic\n stochasticity")),
    list("Deterministic")
  )[[1]]
  
  noise_breaks <- ifelse(
    type_name == "all",
    list(c(-0.75, -0.25, seq(0, 2.0, by = 0.25))),
    list(c(-0.25, seq(0, 2.0, by = 0.25)))
  )[[1]]
  
  value_labels <- ifelse(
    leg_label == "Probability",
    label_percent(),
    function (x) {x}
  )
  
  include_det = type_name == "all"
  
  sigma_width = unique(diff(full_stats_df$sigma))[2]
  R0_height = unique(diff(full_stats_df$R0))[2]
  
  fixed_df = full_stats_df %>%
    filter(name == output_name) %>% 
    dplyr::select(sigma, R0, mean, type) %>% 
    filter(if_all(c(sigma, R0, mean), ~ is.finite(.x) & !is.na(.x)))
  
  # contour_df = smooth_zero_contour(fixed_df, output_name, type_name, zero_proxy)
  smooth_level = 1
  contour_df = fixed_df %>% 
    # Stretch out sigma = 0
    stretch_sigma(., include_det) %>% 
    filter(type == type_name) %>% dplyr::select(-type)
  nn_df = get.knnx(data = cbind(contour_df$sigma, contour_df$R0),
                   query = cbind(contour_df$sigma, contour_df$R0),
                   k = smooth_level)
  
  contour_df$var_smooth <- rowMeans(matrix(contour_df[["mean"]][nn_df$nn.index], ncol = smooth_level))
  
  out = full_stats_df %>% 
    ungroup() %>% 
    # Stretch out sigma = 0
    stretch_sigma(., include_det) %>% 
    filter(name == output_name, type == type_name) %>%
    ggplot(aes(y = R0 - R0_height/2)) +
    # Plot values
    geom_tile(aes(
      x = sigma - sigma_width/2,
      fill = mean,
      width = sigma_width,
      height = R0_height
    )) +
    # Add smoothed contour for zero proxy values
    geom_hline(yintercept = 1, color = "red", lwd = 1) +
    geom_contour(
      data = contour_df,
      aes(x = sigma,
          y = R0,
          z = var_smooth),
      breaks = c(zero_proxy),
      color = "white"
    ) +
    # geom_path(data = contour_df,
    #           aes(
    #             x = sigma,
    #             y = R0,
    #             group = group),
    #           color = "white", lwd = 0.5, alpha = 0.75) +
    geom_vline(xintercept = 0, color = "grey", lwd = 1) +
    # Annotate x-axis for demographic noise
    scale_x_continuous(
      TeX("Environmental noise strength [$\\sigma$]"), 
      # limits = c(-0.25, NA),
      expand = c(0,0),
      breaks = noise_breaks,
      labels = c(noise_labels, seq(0, 2.0, by = 0.25))
    ) +
    scale_y_continuous(
      name = TeX("             Basic reproduction number  [$R_0$]"),
      limits = c(0, 5),
      expand = c(0,0)
    ) +
    # color:
    scale_fill_viridis_c(
      name = output_name,
      option = palette,
      limits = case_when(
        output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ c(0,10),
        output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ c(0,1),
        output_name %in% c("max_cases") ~ c(0,10000)
        # TRUE ~ c(0,NA),
      ),
      breaks = case_when(
        output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ c(0, 5, 10),
        output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ c(0,0.5, 1),
        output_name %in% c("max_cases") ~ c(0,5000, 10000)
        # TRUE ~ c(0,NA),
      ),
      labels = value_labels
    ) +
    # legend:
    guides(
      fill = guide_colourbar(
        title = leg_label,
        position = "top",
        direction = "horizontal",
        title.position = "left",
        title.hjust = 0,
        title.vjust = 0.9,
        barwidth = 8,
        show.limits = TRUE,
        draw.ulim = TRUE,
        draw.llim = TRUE,
      )
    ) +
    # ggtitle(nice_output_labeller(output_name)) +
    theme_half_open(8) +
    theme(
      legend.position = "right",
      legend.justification = "right",
      legend.box.just = "right"
    )
  
  if (type_name == "all") {
    out = out + geom_vline(xintercept = -0.5, color = "grey", lwd = 1)
  }
  
  out
  
  
  return(out)
}


# TODO
# [] Set up diverging palettes with set points at the maximum, zero, and minimum
# [] Do a better job with "dieout" cases
# [] - Restrict to cases where duration < 10 years and duration > 0 years
#    - Plot peak cases only when outbreak dies out
# [] In Julia: define outbreak duration as max(time[H > 1]) - min(time[H > 1])
# [] Make more noise difference comparison plots
#    - All noise vs. no environmental (sigma = 0) NB: this is already available
#    - Enviro only vs. deterministic
#    - Dem only vs. deterministic

# Figure 2: Example trajectories ----

# [] Filter this out more. Either:
# [] 1. choose fewer trajectories
# [] 2. plot rolling averages

# x-axis: Time [years]
# y-axis: Number of hosts infected
# Rows = R0 values
# Columns = sigma values

R0s_to_plot = c(0.95, 1.05, 2, 4)  #c(1.125, 1.375, 3, 4.625) #c(0.95, 1.05, 1.25, 2)
sigmas_to_plot = c(0, 0.5, 1, 1.5) #c(0.65, 1, 1.5, 1.65)#c(0.05, 0.1, 0.25, 0.4)

quadrant_df <- All_sims_plot_df %>%
  filter(R0_factor %in% R0s_to_plot,
         sigma_factor %in% sigmas_to_plot)

deterministic_df <- read_csv("./data/trajectories_for_grid_plot_det.csv.gz") %>% 
  mutate(R0 = R0_from_Thv_function(Thv)) %>% 
  dplyr::select(-V)
deterministic_df$R0_factor = factor(round(deterministic_df$R0,3), levels = rev((unique(round(deterministic_df$R0,3)))))
deterministic_df <- deterministic_df %>% 
  filter(R0_factor %in% R0s_to_plot) %>%
  dplyr::select(-sigma) %>% 
  expand_grid(.,data.frame(sigma = sigmas_to_plot))
deterministic_df$sigma_factor = factor(round(deterministic_df$sigma,3), levels = unique(round(deterministic_df$sigma,3)))

summary_df <- quadrant_df  %>%
  dplyr::select(sigma_factor, R0_factor, run, endemic) %>%
  distinct() %>%
  group_by(sigma_factor, R0_factor) %>%
  summarize(prop_died = 100*(max(quadrant_df$run) - sum(endemic))/max(quadrant_df$run)) %>%
  mutate(label = paste0("Proportion died out = ", prop_died, "%")) %>%
  mutate(x = Inf, y = -Inf)


library(slider)
quadrant_plot <- quadrant_df %>%
  filter(run < 21) %>%
  # Plot 30-day moving averages
  arrange(time) %>% 
  group_by(run, Thv, sigma) %>% 
  mutate(moving_avg = slide_dbl(H, mean, .before= 30, .complete = T)) %>% 
  ggplot() +
  # Plot stochastic trajectories
  geom_line(aes(x = (time/365), y = H, color = endemic,
                group = as.factor(run)),
            alpha = 0.25) +
  # Plot mean-field deterministic trajectories
  geom_line(data = deterministic_df,
            aes(x = time/365, y = H)) +
  # Indicate number that died out
  geom_label(data = summary_df, aes(x = x, y = y, label = label),
             label.size = NA, size = 2, hjust = "inward", vjust = "inward") +
  # Color trajectories according to whether they eventually go to zero (=="plum4")
  scale_color_manual(
    name = "Did the outbreak last 10 years?",
    values = c4a("brewer.dark2", 2),
    breaks = c("TRUE", "FALSE"),
    labels = c("Yes", "No")
  ) +
  # Subplots for each combination of environmental noise strength and R0
  facet_grid(R0_factor ~ sigma_factor,
             labeller = labeller(.rows = as_labeller(appender_R0,
                                                     default = label_parsed),
                                 .cols = as_labeller(appender_sigma,
                                                     default = label_parsed)),
             scales = "free_y") +
  scale_x_continuous("Time [years]",
                     breaks = seq(0, 10),
                     expand = c(0,0)) +
  scale_y_continuous("Number of cases in hosts",
                     expand = c(0,0)) +
  theme_minimal_grid(10) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill=NA)
  ) +
  guides(
    color = guide_legend(override.aes = list(lwd = 2, size = 2, alpha = 1))
  )

quadrant_plot

ggsave("./figures/Figure2.png", quadrant_plot, width = 6.5, height = 5, units = "in", dpi = 1200)


# Figure 3: Mean prob/intensity/duration vs. sigma curves ----



# x: sigma
# y: value

# Rows = type
# Color = R0 value

plot_df <- all_stats_df %>% 
  mutate(round_sigma = round(sigma, 3)) %>%
  filter(
    name %in% c("big_outbreak", "max_cases", "duration"),
    R0_factor %in% c(0.95, 1.05, 1, 2, 3, 4, 5),#c(seq(0.25, 1.5, by = 0.25), 2.0, 3.0, 4.0, 5.0),
    round_sigma %in% round(seq(0.0, 2.0, by = 0.05), 3)
  ) %>% 
  mutate(label = case_when(name == "big_outbreak" ~ "A. Probability",
                           name == "max_cases" ~ "B. Intensity",
                           name == "duration" ~ "C. Duration"))

# R0_colors = c(
#   rev(c4a("kovesi.bu_bk_rd", 10)[5:10]),
#   "black",
#   rev(c4a("kovesi.bu_bk_rd", 10)[1:1])
# )
R0_colors = rev(c4a("kovesi.bu_bk_rd", 7))


Big_outbreak_mean <- plot_df %>% 
  arrange(R0) %>% 
  filter(name == "big_outbreak") %>% 
  ggplot(aes(x = sigma)) +
  geom_line(aes(y = mean, color = R0_factor, group = R0_factor)) +
  geom_point(aes(y = mean, color = R0_factor, group = R0_factor)) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = R0_factor), alpha = 0.1) +
  scale_color_manual(values = R0_colors) + #, breaks = R0_vals) +
  scale_fill_manual(values = R0_colors)  +#, breaks = R0_vals) +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  ggtitle("A. Probability") +
  scale_x_continuous(
    unname(TeX("Environmental noise strength $[\\sigma]$")),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Probability",
    labels = scales::percent_format(accuracy = 1),
    expand = c(0,0)
  ) +
  theme(strip.text.x = element_text(hjust = 0),
        strip.background = element_blank()) +
  theme_minimal_grid(10) +
  guides(
    color = guide_none(), fill = guide_none()
  )

Peak_cases_mean <- plot_df %>% 
  filter(name == "max_cases") %>% 
  ggplot(aes(x = sigma)) +
  geom_line(aes(y =  mean, color = R0_factor, group = R0_factor)) +
  geom_point(aes(y =  mean, color = R0_factor, group = R0_factor)) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = R0_factor), alpha = 0.1) +
  scale_color_manual(values = R0_colors) + #, breaks = R0_vals) +
  scale_fill_manual(values = R0_colors)  +#, breaks = R0_vals) +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  ggtitle("B. Intensity") +
  scale_x_continuous(
    unname(TeX("Environmental noise strength $[\\sigma]$")),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Number of cases",
    breaks = c(2500, 5000, 7500, 10000),
    labels = c(2500, 5000, 7500, 10000),
    expand = c(0,400)
  ) +
  theme(strip.text.x = element_text(hjust = 0),
        strip.background = element_blank()) +
  theme_minimal_grid(10)

Duration_mean <- plot_df %>% 
  filter(name == "duration") %>% 
  ggplot(aes(x = sigma)) +
  geom_line(aes(y =  mean, color = R0_factor, group = R0_factor)) +
  geom_point(aes(y =  mean, color = R0_factor, group = R0_factor)) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = R0_factor), alpha = 0.1) +
  scale_color_manual(values = R0_colors) + #, breaks = R0_vals) +
  scale_fill_manual(values = R0_colors)  +#, breaks = R0_vals) +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  ggtitle("C. Duration") +
  scale_x_continuous(
    unname(TeX("Environmental noise strength $[\\sigma]$")),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Time [years]",
    breaks = c(2.5, 5.0, 7.5, 10.0),
    expand = c(0,0.25)
  ) +
  theme(strip.text.x = element_text(hjust = 0),
        strip.background = element_blank()) +
  theme_minimal_grid(10) +
  guides(
    color = guide_none(), fill = guide_none()
  )

Figure3 = egg::ggarrange(Big_outbreak_mean, Peak_cases_mean, Duration_mean, ncol = 1)
ggsave("./figures/Figure3.png", Figure3, width = 6.5, height = 7, units = "in", dpi = 1200)


# Figure 4: Stacked heatmaps with all noise ----
# Pr(> 100 hosts), peak # of cases, Duration


Big_outbreak_heat <- generic_heat_function("big_outbreak", "all") +
  labs(
    title = "A."
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_blank()
  ) +
  theme(legend.margin=margin(t=-0.5,l=0.0,b=-0.275,r=0.25, unit='cm'))


Peak_cases_heat <- generic_heat_function("max_cases", "all") +
  labs(
    title = "B."
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_text(hjust = 0.9)
  ) +
  theme(legend.margin=margin(t=-0.75,l=0.0,b=-0.275,r=0.25, unit='cm'))

Duration_heat <- generic_heat_function("duration", "all") +
  labs(
    title = "C."
  ) +
  theme(
    axis.title.y = element_blank()
  ) +
  theme(legend.margin=margin(t=-0.75,l=0.0,b=-0.275,r=0.25, unit='cm'))

Figure4 = egg::ggarrange(Big_outbreak_heat, Peak_cases_heat, Duration_heat, ncol = 1)
ggsave("./figures/Figure4.png", Figure4, width = 6.5, height = 4, units = "in", dpi = 1200)
#####

Big_outbreak_heat_alt <- generic_heat_function("big_outbreak", "all") +
  labs(
    title = "A."
  ) +
  theme(
    axis.text.y.right = element_blank(),
  ) +
  theme(legend.margin=margin(t=-0.75,l=0.0,b=-0.25,r=0.25, unit='cm'))


Peak_cases_heat_alt <- generic_heat_function("max_cases", "all") +
  labs(
    title = "B."
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),  # Remove x-axis line
    axis.text.y.right = element_blank(),
  ) +
  theme(legend.margin=margin(t=-0.75,l=0.0,b=-0.25,r=0.25, unit='cm'))

Duration_heat_alt <- generic_heat_function("duration", "all") +
  labs(
    title = "C."
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),  # Remove x-axis line
  ) +
  theme(legend.margin=margin(t=-0.75,l=0.0,b=-0.25,r=0.25, unit='cm'))


Figure4_alt = egg::ggarrange(Big_outbreak_heat_alt, Peak_cases_heat_alt, Duration_heat_alt, nrow = 1)
ggsave("./figures/Figure4_alt.png", Figure4_alt, width = 16, height = 8, units = "in", dpi = 1200)

# Figure 5: Stacked heatmaps with no demographic noise ----
# Pr(> 100 hosts), peak # of cases, Duration
Big_outbreak_heat <- generic_heat_function("big_outbreak", "enviro") +
  labs(
    title = "A."
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_blank()
  ) +
  theme(legend.margin=margin(t=-0.5,l=0.0,b=-0.275,r=0.25, unit='cm'))


Peak_cases_heat <- generic_heat_function("max_cases", "enviro") +
  labs(
    title = "B."
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_text(hjust = 0.9)
  ) +
  theme(legend.margin=margin(t=-0.75,l=0.0,b=-0.275,r=0.25, unit='cm'))

Duration_heat <- generic_heat_function("duration", "enviro") +
  labs(
    title = "C."
  ) +
  theme(
    axis.title.y = element_blank()
  ) +
  theme(legend.margin=margin(t=-0.75,l=0.0,b=-0.275,r=0.25, unit='cm'))

Figure5 = egg::ggarrange(Big_outbreak_heat, Peak_cases_heat, Duration_heat, ncol = 1)
ggsave("./figures/Figure5.png", Figure5, width = 6.5, height = 4, units = "in", dpi = 1200)

# Comparison plots ----

# Make absolute difference comparison heatmap
compare_heat_function <- function(output_name, in_df, type) {
  
  if (type == "abs_diff") {
    leg_label = paste0( "Difference in ", case_when(
      output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ "probability",
      output_name %in% c("duration", "peak_time", "duration_dieout") ~ "time [years]",
      output_name %in% c("max_cases") ~ "number of cases",
    )
    )
    
    type_label = "absolute difference (w vs. w/o dem. noise)"
    
    # Get number of unique values to assign colors to
    z_vals = in_df %>% 
      ungroup() %>% 
      filter(name == output_name) %>% 
      dplyr::select(abs_diff) %>% unique() %>% 
      pull()
    num_cols = length(z_vals)
    z_min = min(z_vals, na.rm = T)
    z_max = max(z_vals, na.rm = T)
    act_min = min(z_min, -z_max)
    act_max = max(z_max, -z_min)
    mid_val = 0
    
    # Set heatmap palette to match unit of output
    palette = case_when(
      output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ c4a("hcl.red_green", num_cols, type = "div"),
      output_name %in% c("duration", "peak_time", "duration_dieout") ~ rev(c4a("matplotlib.seismic", num_cols, type = "div")),
      output_name %in% c("max_cases") ~ c4a("cols4all.pu_gn_div", num_cols, type = "div")	
    )
    # palettes with zero = black
    # kovesi.bu_bk_br
    # scico.managua	
    # scico.vanimo
    
    # palettes with zero = white
    # hcl.red_green matplotlib.seismic cols4all.pu_gn_div
    
  } else {
    leg_label = "% difference"
    in_df <- mutate(in_df, perc_diff = 100 * perc_diff)
    palette = "magma"
    
    type_label = "percent difference (w vs. w/o dem. noise)"
  }
  
  # Approximate difference as a continuous function to get smooth contour lines
  fixed_df = in_df %>%
    filter(type == "all", name == output_name) %>% 
    dplyr::select(sigma, R0, abs_diff) %>% 
    filter(if_all(c(sigma, R0, abs_diff), ~ is.finite(.x) & !is.na(.x))) %>% 
    arrange(sigma, R0)
  
  sigma_width = unique(diff(fixed_df$sigma))[2]
  R0_height = unique(diff(fixed_df$R0))[2]
  
  # contour_df <- smooth_zero_contour(fixed_df, output_name, "comp", 0.01) # (in_df, output_name, type_name, level_val)
  smooth_level = 10
  contour_df = fixed_df %>% 
    stretch_sigma()
  
  nn_df = get.knnx(data = cbind(contour_df$sigma, contour_df$R0),
                   query = cbind(contour_df$sigma, contour_df$R0),
                   k = smooth_level)
  
  contour_df$var_smooth <- rowMeans(matrix(contour_df[["abs_diff"]][nn_df$nn.index], ncol = smooth_level))
  
  
  fixed_df %>% 
    # Stretch out sigma = 0
    stretch_sigma() %>%
    ggplot() +
    # Tile values across grid
    geom_tile(aes(
      x = sigma - sigma_width/2,
      y = R0 - R0_height/2,
      fill = !!sym(type)
    )) +
    # Add a red line for R0 = 1
    geom_hline(yintercept = 1, color = "grey", lwd = 1) +
    # Add a black line for "no environmental noise"
    geom_vline(xintercept = 0, color = "grey", lwd = 1) +
    # Add smoothed contour for abs_diff = 0
    geom_contour(
      data = contour_df,
      aes(
        x = sigma,
        y = R0,
        z = var_smooth),
      breaks = c(0.01),
      color = "orange"
    ) +
    # geom_path(data = contour_df,
    #           aes(
    #             x = sigma - sigma_width,
    #             y = R0 + R0_height,
    #             group = group),
    #           color = "black", lwd = 0.5, alpha = 0.75) +
    # Annotate x-axis for demographic noise
    scale_x_continuous(TeX("Environmental noise strength [$\\sigma$]"),
                       limits = c(-0.5, 2),
                       breaks = c(-0.24, seq(0, 2.0, by = 0.25)),
                       labels = c("Difference from\n deterministic submodel", seq(0, 2.0, by = 0.25)),
                       expand = c(0,0),
    ) +
    scale_y_continuous(TeX("Basic reproduction number  [$R_0$]"),
                       limits = c(0, 5),
                       expand = c(0,0)) +    
    # color:
    scale_fill_gradientn(
      colors = palette,
      limits = c(z_min, z_max),
      values = scales::rescale(c(z_min, mid_val, z_max)),
      labels = ifelse(
        output_name %in% c("small_outbreak", "big_outbreak", "endemic"),
        function(x) paste0(100*x, "%"),
        function(x) x)
      # oob = scales::squish
    ) +
    # legend:
    guides(
      fill = guide_colourbar(
        title = leg_label,
        position = "top",
        direction = "horizontal",
        title.position = "left",
        title.hjust = 0,
        title.vjust = 0.9,
        barwidth = 8,
        show.limits = TRUE,
        draw.ulim = TRUE,
        draw.llim = TRUE,
      )
    ) +
    # ggtitle(nice_output_labeller(output_name)) +
    theme_half_open(8) +
    theme(
      legend.position = "right",
      legend.justification = "right",
      legend.box.just = "right",
      axis.text.y.right = element_text(size = 8)
    )
}

# Figure 6: stacked comparison heatmaps ----
Big_outbreak_plot <- compare_heat_function("big_outbreak", comp_stats_df, "abs_diff") +
  labs(
    title = "A."
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_blank()
  ) +
  theme(legend.margin=margin(t=-0.5,l=0.0,b=-0.35,r=0, unit='cm'))

Peak_cases_plot <- compare_heat_function("max_cases", comp_stats_df, "abs_diff") +
  labs(
    title = "B."
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_text(hjust = 0.9)
  ) + 
  theme(legend.margin=margin(t=-0.75,l=0.0,b=-0.275,r=0, unit='cm'))

Duration_plot <- compare_heat_function("duration", comp_stats_df, "abs_diff") +
  labs(
    title = "C."
  ) +
  theme(
    axis.title.y = element_blank()
  ) +
  theme(legend.margin=margin(t=-0.75,l=0.0,b=-0.275,r=0, unit='cm'))

Figure6 = egg::ggarrange(Big_outbreak_plot, Peak_cases_plot, Duration_plot, ncol = 1)
ggsave("./figures/Figure6.png", Figure6, width = 6.5, height = 4, units = "in", dpi = 1200)

# New Fig 3: Outbreak probability plots ----

# [] Combine line plot and heatmaps (complete model and environmental noise submodel)
# [] Stack heatmaps on right. Put line plot on left matching the height of heatmaps
# [] Use a single color bar for the heatmaps and position it appropriately
# [] Label the subplots A. B., and C. with identifying subtitles
# [] Remove "Deterministic submodel" parts from the heatmaps




# Supp Fig 1: Peak case count explanatory plot ----

# Plot trying to explain difference in peak cases between all noise and no demographic noise
# Find coordinates of the "cool spot" in the peak cases plots
peak_coolspot = comp_stats_df %>% 
  # ungroup() %>% 
  filter(name == "max_cases") %>% 
  filter(abs_diff < 0.99 * min(abs_diff))

peak_point_plot <- comparison_trajectories %>% 
  mutate(
    alpha_level = case_when(
      type %in% c("All_noise", "No_demographic") ~ 0.00001,
      TRUE ~ 1
    )
  ) %>% 
  ggplot(aes(color = type)) +
  # Solution trajectories
  # geom_line(aes(x = time, y = H, group = interaction(R0, sigma, type, run), alpha = alpha_level)) +
  # Mean of stochastic trajectories
  geom_line(aes(x = time, y = mean_H, group = interaction(R0, sigma, type))) +
  # Peak points
  geom_point(
    data = comparison_trajectories %>% 
      dplyr::select(max_time, max_H, type, R0, sigma) %>% 
      unique(),
    aes(x = max_time, y = max_H, group = interaction(type, R0, sigma)),
    alpha = 0.5) +
  # Subplots for each combination of environmental noise strength and R0
  facet_grid(R0 ~ sigma,
             labeller = labeller(.rows = as_labeller(appender_R0,
                                                     default = label_parsed),
                                 .cols = as_labeller(appender_sigma,
                                                     default = label_parsed)),
             scales = "free_y") +
  scale_color_manual(
    values = c4a("met.lakota",3)
  ) + 
  guides(
    alpha = "none"
  ) +
  ggtitle("Points show Intensity for each run of the stochastic models") +
  theme_minimal(11) +
  theme_minimal(11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"    
  )


peak_point_plot
ggsave("./figures/peak_cases_comparison.png", peak_point_plot, width = 6.5, height = 4.5, units = "in")

# Supp Fig 2: Peak histogram plot ----

peak_histogram_plot_no_condition <- comparison_trajectories %>% 
  filter(type != "Deterministic") %>% 
  dplyr::select(max_time, max_H, type, R0, sigma) %>% 
  # filter(max_H > 100) %>% 
  unique() %>% 
  ggplot(aes(x = max_H, color = type, fill = type)) +
  # Peak points histogram
  geom_histogram(
    aes(y = after_stat(count)),  # Ensure y-axis represents counts per facet
    # position = "identity",  # Prevent stacking
    position = 'dodge',
    alpha = 0.5,
    # bins = 30  # Set fixed bin count (adjustable)
  ) +
  # Means of distributions
  geom_vline(
    data = comparison_trajectories %>% 
      filter(type != "Deterministic") %>% 
      dplyr::select(max_time, max_H, type, R0, sigma) %>% 
      # filter(max_H > 100) %>% 
      unique() %>% 
      group_by(type, R0, sigma) %>% 
      summarise(mean_max_H = mean(max_H)),
    aes(xintercept = mean_max_H, color = type),
    lwd = 1, lty = 2
  ) +
  # Subplots for each combination of environmental noise strength and R0
  ggh4x::facet_grid2(R0 ~ sigma,
                     scales = "free",
                     independent = "all",
                     labeller = labeller(.rows = as_labeller(appender_R0,
                                                             default = label_parsed),
                                         .cols = as_labeller(appender_sigma,
                                                             default = label_parsed))
  ) +
  scale_color_manual(
    values = c4a("met.lakota",3)[2:3]
  ) + 
  scale_fill_manual(
    values = c4a("met.lakota",3)[2:3]
  ) + 
  guides(
    alpha = "none"
  ) +
  theme_minimal(11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"    
  )

peak_histogram_plot_no_condition
ggsave("./figures/peak_cases_comparison_histogram_no_condition.png", peak_histogram_plot_no_condition, width = 6.5, height = 4.5, units = "in")

# Small outbreaks case
peak_histogram_plot_small_outbreak <- comparison_trajectories %>% 
  filter(type != "Deterministic") %>% 
  dplyr::select(max_time, max_H, type, R0, sigma) %>% 
  filter(max_H > 10) %>%
  unique() %>% 
  ggplot(aes(x = max_H, color = type, fill = type)) +
  # Peak points histogram
  geom_histogram(
    aes(y = after_stat(count)),  # Ensure y-axis represents counts per facet
    # position = "identity",  # Prevent stacking
    position = 'dodge',
    alpha = 0.5,
    # bins = 30  # Set fixed bin count (adjustable)
  ) +
  # Means of distributions
  geom_vline(
    data = comparison_trajectories %>% 
      filter(type != "Deterministic") %>% 
      dplyr::select(max_time, max_H, type, R0, sigma) %>% 
      filter(max_H > 10) %>%
      unique() %>% 
      group_by(type, R0, sigma) %>% 
      summarise(mean_max_H = mean(max_H)),
    aes(xintercept = mean_max_H, color = type),
    lwd = 1, lty = 2
  ) +
  # Subplots for each combination of environmental noise strength and R0
  ggh4x::facet_grid2(R0 ~ sigma,
                     scales = "free",
                     independent = "all",
                     labeller = labeller(.rows = as_labeller(appender_R0,
                                                             default = label_parsed),
                                         .cols = as_labeller(appender_sigma,
                                                             default = label_parsed))
  ) +
  scale_color_manual(
    values = c4a("met.lakota",3)[2:3]
  ) + 
  scale_fill_manual(
    values = c4a("met.lakota",3)[2:3]
  ) + 
  guides(
    alpha = "none"
  ) +
  theme_minimal(11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"    
  )

peak_histogram_plot_small_outbreak
ggsave("./figures/peak_cases_comparison_histogram_small_outbreak.png", peak_histogram_plot_small_outbreak, width = 6.5, height = 4.5, units = "in")

# Big outbreak case
peak_histogram_plot_big_outbreak <- comparison_trajectories %>% 
  filter(type != "Deterministic") %>% 
  dplyr::select(max_time, max_H, type, R0, sigma) %>% 
  filter(max_H > 100) %>%
  unique() %>% 
  ggplot(aes(x = max_H, color = type, fill = type)) +
  # Peak points histogram
  geom_histogram(
    aes(y = after_stat(count)),  # Ensure y-axis represents counts per facet
    # position = "identity",
    position = 'dodge',  # Prevent stacking
    alpha = 0.5,
    # bins = 30  # Set fixed bin count (adjustable)
  ) +
  # Means of distributions
  geom_vline(
    data = comparison_trajectories %>% 
      filter(type != "Deterministic") %>% 
      dplyr::select(max_time, max_H, type, R0, sigma) %>% 
      filter(max_H > 100) %>%
      unique() %>% 
      group_by(type, R0, sigma) %>% 
      summarise(mean_max_H = mean(max_H)),
    aes(xintercept = mean_max_H, color = type),
    lwd = 1, lty = 2
  ) +
  # Subplots for each combination of environmental noise strength and R0
  ggh4x::facet_grid2(R0 ~ sigma,
                     scales = "free",
                     independent = "all",
                     labeller = labeller(.rows = as_labeller(appender_R0,
                                                             default = label_parsed),
                                         .cols = as_labeller(appender_sigma,
                                                             default = label_parsed))
  ) +
  scale_color_manual(
    values = c4a("met.lakota",3)[2:3]
  ) + 
  scale_fill_manual(
    values = c4a("met.lakota",3)[2:3]
  ) + 
  guides(
    alpha = "none"
  ) +
  theme_minimal(11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"    
  )

peak_histogram_plot_big_outbreak
ggsave("./figures/peak_cases_comparison_histogram_big_outbreak.png", peak_histogram_plot_big_outbreak, width = 6.5, height = 4.5, units = "in")


# Compare how the distribution changes as we filter out smaller outbreaks
peak_comparison_df <- comparison_trajectories %>% 
  filter(type != "Deterministic") %>% 
  filter(R0 == 1.375, sigma == 0.55) %>% 
  mutate(
    max_H_10 = if_else(max_H > 10, max_H, NA),
    max_H_100 = if_else(max_H > 100, max_H, NA)
  ) %>% 
  pivot_longer(cols = c(max_H, max_H_10, max_H_100)) %>% 
  dplyr::select(max_time, name, value, type, R0, sigma) %>% 
  unique()

peak_histogram_plot_comparison <- peak_comparison_df %>% 
  ggplot(aes(x = value, color = type, fill = type)) +
  # Peak points histogram
  geom_histogram(
    aes(y = after_stat(count)),  # Ensure y-axis represents counts per facet
    position = "identity",
    # position = 'dodge',  # Prevent stacking
    alpha = 0.5,
    # bins = 30  # Set fixed bin count (adjustable)
  ) +
  # Means of distributions
  geom_vline(
    data = peak_comparison_df %>% 
      group_by(type, name) %>% 
      summarise(mean_max_H = mean(value, na.rm = T)),
    aes(xintercept = mean_max_H, color = type),
    lwd = 1, lty = 2,
    show.legend = F
  ) +
  # Subplots for each combination of environmental noise strength and R0
  facet_wrap(
    ~name,
    ncol = 1,
    # scales = "free",
    # independent = "all",
    labeller = labeller(name = c(
      max_H = "All simulations",
      max_H_10 = "Small outbreaks (>10)",
      max_H_100 = "Large outbreaks (>100)"
    ))
  ) +
  scale_x_continuous(
    name = "Peak case count",
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Count",
    expand = c(0,0)
  ) +
  scale_color_manual(
    name = "Noise type:",
    values = c4a("met.lakota",3)[2:3]
  ) + 
  scale_fill_manual(
    name = "Noise type:",
    values = c4a("met.lakota",3)[2:3]
  ) + 
  guides(
    alpha = "none"
  ) +
  theme_minimal(11) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal"    
  )

peak_histogram_plot_comparison
ggsave("./figures/peak_histogram_plot_comparison.png", peak_histogram_plot_comparison, width = 6.5, height = 4.5, units = "in")


# Figure S3: Scatterplot of dieout duration vs. peak case count ----


outbreak_duration_df <- read_rds("./data/peak_v_duration_sims.rds")

# R0 vals = 0.95, 1.125, 2, 4.625
# sigma vals = 0.25, 0.65, 1, 1.5

duration_peak_scatter_all <- outbreak_duration_df %>% 
  filter(type == "All noise") %>% 
  ggplot(aes(x = max_time, y = max_value)) +
  # geom_hex() +
  geom_pointdensity() +
  geom_point(
    data = outbreak_duration_df %>% group_by(R0_factor, sigma) %>% mutate(mean_max_time = mean(max_time), mean_max_value = mean(max_value)),
    aes(x = mean_max_time, y = mean_max_value),
    color = "red",
    fill = "red",
    shape = 23,
    size = 2
  ) +
  # Subplots for each combination of environmental noise strength and R0
  ggh4x::facet_grid2(
    R0_factor ~ sigma,
    labeller = labeller(.rows = as_labeller(appender_R0, default = label_parsed),
                        .cols = as_labeller(appender_sigma, default = label_parsed))
  ) +
  scale_color_viridis_c(
    name = "Number of\n neighbors"
  ) +
  scale_x_continuous(
    "Duration [years]",
    expand = c(0,0.25),
  ) +
  scale_y_continuous(
    "Intensity [cases]",
    limits = c(0,10000)
  ) +
  # legend:
  theme_minimal(11) +
  theme(
    strip.background = element_rect(color = "white", fill = "white")
  )

duration_peak_scatter_all
ggsave("./figures/duration_peak_scatter_all.png", duration_peak_scatter_all, width = 6.5, height = 4.5, units = "in")

duration_peak_scatter_enviro <- outbreak_duration_df %>% 
  filter(type == "Environmental noise only") %>% 
  ggplot(aes(x = max_time, y = max_value)) +
  # geom_hex() +
  geom_pointdensity() +
  # Subplots for each combination of environmental noise strength and R0
  ggh4x::facet_grid2(
    R0_factor ~ sigma,
    labeller = labeller(.rows = as_labeller(appender_R0, default = label_parsed),
                        .cols = as_labeller(appender_sigma, default = label_parsed))
  ) +
  scale_color_viridis_c(
    name = "Number of\n neighbors"
  ) +
  scale_x_continuous(
    "Duration [years]",
    expand = c(0,0.25),
  ) +
  scale_y_continuous(
    "Intensity [cases]",
    limits = c(0,10000)
  ) +
  # legend:
  theme_half_open(11) +
  theme(
    strip.background = element_rect(color = "white", fill = "white")
  )

duration_peak_scatter_enviro
ggsave("./figures/duration_peak_scatter_enviro.png", duration_peak_scatter_enviro, width = 6.5, height = 4.5, units = "in")


# Figure SXXX: "Bifurcation" diagram of R0 vs. peak cases

# "Bifurcation" diagram of peak cases
bif_plot <- outbreak_duration_df %>% 
  filter(
    type == "Environmental noise only",
    round(sigma, 3) %in% c(0.25, 0.75, 1.25, 1.75),
  ) %>% 
  arrange(sigma) %>% 
  group_by(R0, sigma) %>% 
  mutate(
    mean_max_cases = mean(max_value),
    # dens = approxfun(density(max_cases, kernel = "cosine"))(max_cases)
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = R0, group = as.factor(sigma))) +
  geom_pointdensity(aes(y = max_value)) +
  geom_path(aes(y = mean_max_cases), color = "red", lwd = 1) +
  facet_wrap( ~ as.factor(sigma),
              labeller = as_labeller(appender_sigma, default = label_parsed),
              ncol = 1, scales = "free") +
  scale_color_viridis_c(
    "Number of\nneighboring\npoints"
  ) +
  scale_x_continuous(
    "Basic reproduction number [R0]",
  ) +
  scale_y_continuous(
    "Intensity [cases]",
    expand = c(0,250)
  ) +
  theme_minimal(11) +
  theme(
    strip.background = element_rect(color = "white", fill = "white")
  )

ggsave("./figures/peak_cases_density_bifurcation_plot_switch.png", bif_plot,
       width = 13.5, height = 7.5, units = "in")


# For each of these points, plot distributions of outputs, see how well a Gamma distribution fits.

# Supp Fig XXX: Duration of dieouts Histogram ----
duration_dieout_histograms <- enviro_outbreak_duration_df %>% 
  # filter(!is.na(duration_dieout)) %>% 
  ggplot(aes(x = duration_dieout, group = interaction(R0_factor, sigma))) +
  geom_histogram() +
  # Subplots for each combination of environmental noise strength and R0
  ggh4x::facet_grid2(R0_factor ~ sigma,
                     labeller = labeller(.rows = as_labeller(appender_R0,
                                                             default = label_parsed),
                                         .cols = as_labeller(appender_sigma,
                                                             default = label_parsed)),
                     scales = "free",
                     independent = "all") +
  scale_color_manual(
    values = c4a("met.lakota",3)[2:3]
  ) +
  theme_minimal()



# Fit Gamma distribution, perform tests, and generate plot

fitted_data <- outbreak_duration_df %>%
  group_by(type) %>% 
  filter(!is.na(duration_dieout)) %>% 
  filter(duration_dieout > 1/365) %>%
  dplyr::select(R0_factor, sigma, duration_dieout) %>% 
  rename(value = duration_dieout) %>% 
  group_by(R0_factor, sigma) %>%
  mutate(
    fit = list(tryCatch({
      fit <- fitdistr(value, "gamma")
      list(
        shape = fit$estimate["shape"],
        rate = fit$estimate["rate"],
        p_value = ks.test(value, "pgamma", shape = fit$estimate["shape"], rate = fit$estimate["rate"])$p.value
      )
    }, error = function(e) NULL))
  ) %>%
  unnest_wider(fit) %>%
  filter(!is.na(shape)) # Remove combinations where fitting failed


# Add theoretical density values for Gamma distribution
density_data <- fitted_data %>%
  group_by(R0_factor, sigma) %>%
  summarize(
    density_x = seq(min(value), max(value), length.out = 100),
    density_y = dgamma(seq(min(value), max(value), length.out = 100), shape = unique(shape), rate = unique(rate)),
    .groups = "drop"
  )

# Create ggplot with faceting
fitted_data %>% 
  # filter(value > 1 / 365) %>% 
  ggplot(aes(x = value)) +
  # Histogram
  geom_histogram(
    aes(y = ..count.., fill = type), 
    position = 'identity',
    # bins = 30, 
    # fill = "blue", 
    alpha = 0.6, 
    # color = "black",
    color = NA
  ) +
  # Means of distributions
  geom_vline(
    data = fitted_data %>%
      group_by(R0_factor, sigma,type) %>%
      summarise(mean_value = mean(value)),
    aes(xintercept = mean_value, color = type),
    lwd = 1, lty = 2
  ) +
  # Subplots for each combination of environmental noise strength and R0
  ggh4x::facet_grid2(R0_factor ~ sigma,
                     labeller = labeller(.rows = as_labeller(appender_R0,
                                                             default = label_parsed),
                                         .cols = as_labeller(appender_sigma,
                                                             default = label_parsed)),
                     scales = "free_y",
                     independent = "y"
  ) +
  scale_x_continuous(
    name = "Outbreak duration [years]",
    trans = 'log10',
    expand = c(0,0)
    # breaks = seq(0, 10)
  ) +
  scale_y_continuous(
    name = "Count",
    expand = c(0,0)
    # trans = 'log10'
    # limits = c(0,1)
  ) +
  labs(
    title = "Histograms with Fitted Gamma Distributions",
    x = "Data Values",
    y = "Density"
  ) +
  theme_minimal(11) +
  # theme(strip.text = element_text(size = 11), panel.spacing = unit(1, "lines")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"    
  )


ggsave("./figures/duration_dieout_fits.png", width = 6.5, height = 4.5, units = "in")


# Set heatmap palette to match unit of output
palette = "plasma"

# Label legend appropriately
leg_label = "Time [years]"

# noise_label = case_when(
#   type_name == "all" ~  "Demo\n noise only",
#   type_name == "enviro" ~ "No\n noise"
# )

# 
# all_df_modified %>% 
#   ungroup() %>% 
#   # Stretch out sigma = 0
#   filter(
#     !is.na(mean),
#     name == "duration_10"
#   ) %>%
#   # stretch_sigma() %>% 
#   ggplot(aes(x = sigma, y = R0, z = mean)) +
#   geom_raster(aes(fill = mean),
#               hjust = unique(diff(all_df_modified$sigma))[2] / 2,
#               vjust = unique(diff(all_df_modified$R0))[2] / 2) +
#   geom_hline(yintercept = 1, color = "red", lwd = 1) +
#   geom_vline(xintercept = 0, color = "black", lwd = 1) +
#   # Annotate x-axis for demographic noise
#   scale_x_continuous(TeX("Environmental noise strength [$\\sigma$]"),
#                      limits = c(-0.25, NA),
#                      expand = c(0,0),
#                      breaks = c(-0.125, seq(0, 2.0, by = 0.25)),
#                      labels = c("Demo\n noise only", seq(0, 2.0, by = 0.25))
#   ) +
#   scale_y_continuous("R0",
#                      limits = c(0, NA),
#                      expand = c(0,0)) +    
#   # color:
#   scale_fill_viridis_c(
#     name = "Outbreak duration (not endemic, > 10 cases)",
#     option = "plasma"
#   ) +
#   # legend:
#   guides(fill = guide_colourbar(
#     title = leg_label,
#     title.position = "top",
#     title.hjust = 0.5,
#     barheight = 10,
#     show.limits = TRUE
#   )) +
#   theme_half_open(10)
# 
# outbreak_comparison_df %>% 
#   filter(type == "All" ) %>% 
#   ggplot(aes(x = sigma, y = R0, fill = max_time_10)) +
#   geom_tile()