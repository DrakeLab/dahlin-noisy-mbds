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

# Functions ----

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
    
    # Repeat the zero_sigma rows for each negative sigma value
    det_stretched <- det_zeros[rep(1:nrow(det_zeros), each = length(det_sigmas)), ]
    
    # Assign the new sigma values
    det_stretched$sigma <- rep(det_sigmas, times = nrow(det_zeros))
    det_stretched = mutate(det_stretched, type = "all")
    
    out <- bind_rows(out, det_stretched)
  }
  return(out)
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
    output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ "Time\n[years]",
    output_name %in% c("max_cases") ~ "Number\nof cases",
  )
  
  # Define zero proxy values
  zero_proxy = case_when(
    output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ 1E-4,
    output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ 1,
    output_name %in% c("max_cases") ~ 1,
  )
  
  noise_labels = ifelse(
    type_name == "all",
    "Demo. noise\nsubmodel",
    "Deterministic\nsubmodel"
  )[[1]]
  
  noise_breaks <- c(-0.25, seq(0, 2.0, by = 0.25))
  
  value_labels <- ifelse(
    leg_label == "Probability",
    label_percent(),
    function (x) {x}
  )
  
  sigma_width = unique(diff(full_stats_df$sigma))[2]
  R0_height = unique(diff(full_stats_df$R0))[2]
  
  fixed_df = full_stats_df %>%
    filter(name == output_name) %>% 
    dplyr::select(sigma, R0, mean, type) %>% 
    filter(if_all(c(sigma, R0, mean), ~ is.finite(.x) & !is.na(.x)))
  
  smooth_level = 10
  contour_df = fixed_df %>% 
    stretch_sigma()
  
  nn_df = get.knnx(data = cbind(contour_df$sigma, contour_df$R0),
                   query = cbind(contour_df$sigma, contour_df$R0),
                   k = smooth_level)
  
  contour_df$var_smooth <- rowMeans(matrix(contour_df[["mean"]][nn_df$nn.index], ncol = smooth_level))
  
  out = full_stats_df %>% 
    ungroup() %>% 
    filter(name == output_name, type == type_name) %>%
    stretch_sigma() %>% 
    ggplot() +
    # Plot values
    geom_tile(aes(
      x = sigma - sigma_width/2,
      y = R0 - R0_height/2,
      fill = mean,
      width = sigma_width,
      height = R0_height
    )) +
    # Add smoothed contour for zero proxy values
    geom_hline(yintercept = 1, color = "red", lwd = 0.5) +
    geom_contour(
      aes(
        x = sigma,
        y = R0,
        z = mean),
      breaks = c(zero_proxy),
      color = "grey100"
    ) +
    # Line showing divide from submodel
    geom_vline(
      xintercept = 0,
      color = "grey",
      lwd = 0.5
    ) +
    # Annotate x-axis for demographic noise
    scale_x_continuous(
      TeX("Environmental noise [$\\sigma$]"), 
      expand = c(0,0),
      breaks = noise_breaks,
      labels = c(noise_labels, seq(0, 2.0, by = 0.25))
    ) +
    scale_y_continuous(
      name = ifelse(type_name == "enviro",
                    TeX("                                            Basic reproduction number  [$R_0$]"),
                    "          "),
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
      ),
      breaks = case_when(
        output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ c(0, 5, 10),
        output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ c(0,0.5, 1),
        output_name %in% c("max_cases") ~ c(0,5000, 10000)
      ),
      labels = value_labels
    ) + 
    coord_cartesian(xlim = c(-0.49, NA)) +
    # legend:
    guides(
      fill = guide_colourbar(
        title = leg_label,
        position = "right",
        direction = "vertical",
        title.position = "top",
        title.hjust = 0,
        title.vjust = 1,
        barheight = 8,
        show.limits = TRUE,
        draw.ulim = TRUE,
        draw.llim = TRUE,
      )
    ) +
    theme_half_open(10) +
    theme(
      legend.position = "right",
      legend.justification = "center",
      legend.box.just = "right"
    )
  
  return(out)
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

# Load data from summarize_solutions.R ----

all_stats_df = read_rds("./data/all_stats.rds")
enviro_stats_df = read_rds("./data/enviro_stats.rds")
full_stats_df <- rbind(
  mutate(all_stats_df, type = "all"),
  mutate(enviro_stats_df, type = "enviro")
)
comp_stats_df = read_rds("./data/comp_stats.rds")
All_sims_plot_df = read_rds("./data/all_sims.rds")
comparison_trajectories = read_rds("./data/comp_trajectories.rds")


R0_colors = rev(c4a("kovesi.bu_bk_rd", 7))
line_plots_df <- all_stats_df %>% 
  mutate(round_sigma = round(sigma, 3)) %>%
  filter(
    name %in% c("big_outbreak", "max_cases", "duration"),
    R0_factor %in% c(0.95, 1.05, 1, 2, 3, 4, 5),
    round_sigma %in% round(seq(0.0, 2.0, by = 0.05), 3)
  ) %>% 
  mutate(label = case_when(name == "big_outbreak" ~ "A. Probability",
                           name == "max_cases" ~ "B. Intensity",
                           name == "duration" ~ "C. Duration"))

# Figure 2: Example trajectories ----

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
  # Subplots for each combination of environmental noise and R0
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


# Figure 3: Probability plots ----

## Subplot A) Line plot
#  = full-size, on left
Big_outbreak_mean <- line_plots_df %>% 
  arrange(R0) %>% 
  filter(name == "big_outbreak") %>% 
  mutate(sep_one = round(R0,2) == 1) %>% 
  ggplot(aes(x = sigma)) +
  geom_line(aes(y = mean, color = R0_factor, group = R0_factor, linetype = factor(sep_one)),
            lwd = 1) +
  scale_linetype_manual( 
    values = c(1, 3)
  ) +
  scale_color_manual(
    values = R0_colors,
    breaks = rev(unique(line_plots_df$R0_factor))
  ) +
  scale_fill_manual(
    values = R0_colors,
    breaks = rev(unique(line_plots_df$R0_factor))
  )  +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  ggtitle("A.") +
  scale_x_continuous(
    unname(TeX("Environmental noise $[\\sigma]$")),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Probability",
    labels = scales::percent_format(accuracy = 1),
    expand = c(0,0.005)
  ) +
  theme_minimal_grid(10) +
  guides(
    color = guide_legend(
      title = TeX("Basic reproduction number $[R_0]$:"),
      title.position = "top",
      label.position = "bottom",
      reverse = TRUE,
      direction = "horizontal",
      nrow = 1,
      override.aes = list(linewidth = 3)
    ),
    linetype = guide_none()
  ) +
  theme(
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 7, margin = margin(t = -1)),
    legend.spacing = unit(0, "pt"),
    legend.title = element_text(size = 8, margin = margin(b = -2)),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",  
    legend.spacing.y = unit(0, "pt"),
    legend.key.width = unit(14, "pt"),
    plot.title.position = "panel",
    legend.box.margin = margin(t = -18, b = -10)
  )

## Subplot B) Heatmap for full model
#  = half-size, on top-right
#  * label clearly with subtitle
#  * center colorbar with below figure
Big_outbreak_heat_temp <- generic_heat_function("big_outbreak", "all")

legend_temp = get_legend(Big_outbreak_heat_temp)

Big_outbreak_heat = Big_outbreak_heat_temp +
  labs(
    title = "B. Full model"
  ) +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_blank()
  ) +
  guides(fill = "none")

## Subplot C) Heatmap for environmental noise sub-model
#  = half-size, on bottom-right
#  * label clearly with subtitle
Big_outbreak_heat_enviro <- generic_heat_function("big_outbreak", "enviro") +
  labs(
    title = "C. Environmental noise submodel"
  ) +
  theme(
    plot.title = element_text(size = 10),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) +
  guides(fill = "none")

Figure3 = ggpubr::ggarrange(Big_outbreak_mean,
                            egg::ggarrange(
                              Big_outbreak_heat, 
                              Big_outbreak_heat_enviro, 
                              ncol = 1,
                              heights = c(0.5,0.5)
                            ),
                            legend_temp,
                            ncol = 3,
                            widths = c(0.4, 0.525, 0.075)
)

ggsave("./figures/Figure3.png", Figure3, width = 9, height = 4, units = "in", dpi = 1200)

# Figure 4: Intensity plots ----

## Subplot A) Line plot
#  = full-size, on left
Peak_cases_mean <- line_plots_df %>% 
  arrange(R0) %>% 
  filter(name == "max_cases") %>% 
  mutate(sep_one = round(R0,2) == 1) %>% 
  ggplot(aes(x = sigma)) +
  geom_line(aes(y = mean, color = R0_factor, group = R0_factor, linetype = factor(sep_one)),
            lwd = 1) +
  scale_color_manual(
    values = R0_colors,
    breaks = rev(unique(line_plots_df$R0_factor))
  ) +
  scale_fill_manual(
    values = R0_colors,
    breaks = rev(unique(line_plots_df$R0_factor))
  )  +
  scale_linetype_manual( 
    values = c(1, 3)
  ) +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  ggtitle("A.") +
  scale_x_continuous(
    unname(TeX("Environmental noise $[\\sigma]$")),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Number of cases",
    breaks = c(2500, 5000, 7500, 10000),
    labels = c(2500, 5000, 7500, 10000),
    expand = expansion(add = c(10,400))
  ) +
  theme_minimal_grid(10) +
  guides(
    color = guide_legend(
      title = TeX("Basic reproduction number $[R_0]$:"),
      title.position = "top",
      label.position = "bottom",
      reverse = TRUE,
      direction = "horizontal",
      nrow = 1,
      override.aes = list(linewidth = 3)
    ),
    linetype = guide_none()
  ) +
  theme(
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 7, margin = margin(t = -1)),
    legend.spacing = unit(0, "pt"),
    legend.title = element_text(size = 8, margin = margin(b = -2)),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",  
    legend.spacing.y = unit(0, "pt"),
    legend.key.width = unit(14, "pt"),
    plot.title.position = "panel",
    legend.box.margin = margin(t = -18, b = -10)
  )

## Subplot B) Heatmap for full model
#  = half-size, on top-right
#  * label clearly with subtitle
#  * center colorbar with below figure
Peak_cases_heat_temp <- generic_heat_function("max_cases", "all")

legend_temp = get_legend(Peak_cases_heat_temp)

Peak_cases_heat = Peak_cases_heat_temp +
  labs(
    title = "B. Full model"
  ) +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_blank()
  ) +
  guides(fill = "none")

## Subplot C) Heatmap for environmental noise sub-model
#  = half-size, on bottom-right
#  * label clearly with subtitle
Peak_cases_heat_enviro <- generic_heat_function("max_cases", "enviro") +
  labs(
    title = "C. Environmental noise submodel"
  ) +
  theme(
    plot.title = element_text(size = 10),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) +
  guides(fill = "none")

Figure4 = ggpubr::ggarrange(Peak_cases_mean,
                            egg::ggarrange(
                              Peak_cases_heat, 
                              Peak_cases_heat_enviro, 
                              ncol = 1,
                              heights = c(0.5,0.5)
                            ),
                            legend_temp,
                            ncol = 3,
                            widths = c(0.4, 0.525, 0.075)
)

ggsave("./figures/Figure4.png", Figure4, width = 9, height = 4, units = "in", dpi = 1200)

# Figure 5: Duration plots ----

## Subplot A) Line plot
#  = full-size, on left
Duration_mean <- line_plots_df %>% 
  arrange(R0) %>% 
  filter(name == "duration") %>% 
  mutate(sep_one = round(R0,2) == 1) %>% 
  ggplot(aes(x = sigma)) +
  geom_line(aes(y = mean, color = R0_factor, group = R0_factor, linetype = factor(sep_one)),
            lwd = 1) +
  scale_color_manual(
    values = R0_colors,
    breaks = rev(unique(line_plots_df$R0_factor))
    ) +
  scale_fill_manual(
    values = R0_colors,
    breaks = rev(unique(line_plots_df$R0_factor))
    )  +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  ggtitle("A.") +
  scale_x_continuous(
    unname(TeX("Environmental noise $[\\sigma]$")),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Time [years]",
    breaks = c(2.5, 5.0, 7.5, 10.0),
    expand = expansion(add = c(0,0.25))
  ) +
  scale_linetype_manual( 
    values = c(1, 3)
    ) +
  theme_minimal_grid(10) +
  guides(
    color = guide_legend(
      title = TeX("Basic reproduction number $[R_0]$:"),
      title.position = "top",
      label.position = "bottom",
      reverse = TRUE,
      direction = "horizontal",
      nrow = 1,
      override.aes = list(linewidth = 3)
    ),
    linetype = guide_none()
  ) +
  theme(
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 7, margin = margin(t = -1)),
    legend.spacing = unit(0, "pt"),
    legend.title = element_text(size = 8, margin = margin(b = -2)),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",  
    legend.spacing.y = unit(0, "pt"),
    legend.key.width = unit(14, "pt"),
    plot.title.position = "panel",
    legend.box.margin = margin(t = -18, b = -10)
  )

## Subplot B) Heatmap for full model
#  = half-size, on top-right
#  * label clearly with subtitle
#  * center colorbar with below figure
Duration_heat_temp <- generic_heat_function("duration", "all")

legend_temp = get_legend(Duration_heat_temp)

Duration_heat = Duration_heat_temp +
  labs(
    title = "B. Full model"
  ) +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_blank()
  ) +
  guides(fill = "none")

## Subplot C) Heatmap for environmental noise sub-model
#  = half-size, on bottom-right
#  * label clearly with subtitle
Duration_heat_enviro <- generic_heat_function("duration", "enviro") +
  labs(
    title = "C. Environmental noise submodel"
  ) +
  theme(
    plot.title = element_text(size = 10),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) +
  guides(fill = "none")

Figure5 = ggpubr::ggarrange(Duration_mean,
                            egg::ggarrange(
                              Duration_heat, 
                              Duration_heat_enviro, 
                              ncol = 1,
                              heights = c(0.5,0.5)
                            ),
                            legend_temp,
                            ncol = 3,
                            widths = c(0.4, 0.525, 0.075)
)

ggsave("./figures/Figure5.png", Figure5, width = 9, height = 4, units = "in", dpi = 1200)

# Figure 6: Comparison plots ------------------------------------------

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
  
  smooth_level = 100
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
      breaks = c(1E-5),
      color = "black"
    ) +
    # Annotate x-axis for demographic noise
    scale_x_continuous(TeX("Environmental noise [$\\sigma$]"),
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
    theme_half_open(8) +
    theme(
      legend.position = "right",
      legend.justification = "right",
      legend.box.just = "right",
      axis.text.y.right = element_text(size = 8)
    )
}

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
ggsave("./figures/Figure6.png", Figure6, width = 4.5, height = 3, units = "in", dpi = 1200)


# Figure S1: Duration - Intensity scatterplots ----------------------------

outbreak_duration_df <- read_rds("./data/peak_v_duration_sims.rds")

# R0 vals = 0.95, 1.125, 2, 4.625
# sigma vals = 0.25, 0.65, 1, 1.5

duration_peak_scatter_all <- outbreak_duration_df %>% 
  filter(type == "Full model") %>% 
  ggplot(aes(x = max_time, y = max_value)) +
  geom_pointdensity() +
  geom_point(
    data = outbreak_duration_df %>% group_by(R0_factor, sigma) %>% mutate(mean_max_time = mean(max_time), mean_max_value = mean(max_value)),
    aes(x = mean_max_time, y = mean_max_value),
    color = "red",
    fill = "red",
    shape = 23,
    size = 2
  ) +
  # Subplots for each combination of environmental noise and R0
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
  ) +
  # legend:
  guides(
    fill = guide_colourbar(
      position = "right",
      direction = "vertical",
      title.position = "top",
      title.hjust = 0,
      title.vjust = 1,
      barheight = 8,
      show.limits = TRUE,
      draw.ulim = TRUE,
      draw.llim = TRUE,
    )
  )

duration_peak_scatter_all
ggsave("./figures/SuppFigure1.png", duration_peak_scatter_all, width = 9, height = 4, units = "in")
