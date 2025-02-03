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


appender_R0 <- function(string) TeX(paste("$R_0 = $", string))
appender_sigma <- function(string) TeX(paste("$\\sigma = $", string))

# # Parameters ---- 
b <- 0.3 # biting rate
Tvh <- 0.5 # vector transmission-probability susceptible vector being infected after contact with infectious host
Nh <- 10000 # number of hosts in system
Nv <- 100000 # number of vectors in system
muv <- 0.1 # vector mortality rate
gammah <- 0.1 # host recovery rate

R0_from_Thv_function <- function(Thv) {
  # Thv <- (R0^2) / ((b^2 * Tvh * Nv) / (Nh * gammah * muv)) 
  sqrt(Thv * b^2 * Tvh * Nv / (Nh * gammah * muv))
}

# Sigma values (environmental noise level in 0-0.3) to test
sigmas <- seq(0, 1, by = 0.05)

# Final time point
max_time = 3650

# Plots ----

# Load data from Julia output
# Deterministic
det_df = read.csv("./data/collect_all_outputs_det.csv") 
# Environmental noise only
enviro_df = read.csv("./data/collect_all_outputs_no_demo.csv") 
# Demographic + Environmental noise
all_df = read.csv("./data/collect_all_outputs.csv") 

# Summarize key statistics: calculate mean, variance, and number of positives
all_summary_df = all_df %>%
  rename(
    max_cases = max_value,
    duration = positive_duration
  ) %>% 
  mutate(
    small_outbreak = ifelse(exceeded_10 == "true", T, F),
    big_outbreak = ifelse(exceeded_100 == "true", T, F),
    endemic = ifelse(positive_at_final == "true", T, F),
    peak_time = max_time / 365, # convert to years
    duration = duration / 365, # convert to years
    duration_dieout = if_else(endemic, NA, 
                              if_else(max_cases <1,
                                      NA,
                                      duration)),
    R0 = R0_from_Thv_function(Thv)
  ) %>% 
  dplyr::select(-c(exceeded_10, exceeded_100, positive_at_final, max_time, Thv)) %>% 
  # Summarise for each combination of R0 and sigma
  pivot_longer(cols = max_cases:duration_dieout) %>%
  mutate(
    .by = c(R0, sigma, name),
    mean = mean(value, na.rm = TRUE),
    variance = var(value, na.rm = TRUE),
    sum_value = as.integer(sum(value, na.rm = TRUE)) # needed for Clopper-Pearon confidence intervals
  ) %>%
  dplyr::select(-c(run, value)) %>% distinct()

enviro_summary_df = enviro_df %>%
  rename(
    max_cases = max_value,
    duration = positive_duration
  ) %>% 
  mutate(
    small_outbreak = ifelse(exceeded_10 == "true", T, F),
    big_outbreak = ifelse(exceeded_100 == "true", T, F),
    endemic = ifelse(positive_at_final == "true", T, F),
    peak_time = max_time / 365, # convert to years
    duration = duration / 365, # convert to years
    duration_dieout = if_else(endemic, NA, duration),
    R0 = R0_from_Thv_function(Thv),
    .keep = "unused"
  ) %>% 
  # Summarise for each combination of R0 and sigma
  pivot_longer(cols = max_cases:duration_dieout) %>%
  mutate(
    .by = c(R0, sigma, name),
    mean = mean(value, na.rm = TRUE),
    variance = var(value, na.rm = TRUE),
    sum_value = as.integer(sum(value, na.rm = TRUE)) # needed for Clopper-Pearon confidence intervals
  ) %>%
  dplyr::select(-c(run, value)) %>% distinct()

num_sims = max(enviro_df$run)

# Process statistics ----

# For both noise types
all_stats_df_binoms <- all_summary_df %>%
  filter(!(name %in% c("max_cases", "duration", "duration_dieout", "peak_time"))) %>%
  rowwise() %>%
  mutate(sum_value = min(sum_value, num_sims)) %>%
  mutate(
    lower_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "lower")$Lower.limit,
    upper_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "upper")$Upper.limit
  ) %>%
  ungroup()

all_stats_df_conts <- all_summary_df %>%
  filter(name %in% c("max_cases", "duration", "duration_dieout", "peak_time")) %>%
  rowwise() %>%
  mutate(
    max_val = if_else(name == "max_cases", Nh, max_time / 365),
    lower_ci = max(0, mean - 1.96 * sqrt(variance)),
    upper_ci = min(max_val, mean + 1.96 * sqrt(variance))
  ) %>%
  dplyr::select(-max_val) %>%
  ungroup()

all_stats_df <- rbind(all_stats_df_binoms, all_stats_df_conts)
all_stats_df$R0_factor = factor(round(all_stats_df$R0,2), levels = rev(unique(round(all_stats_df$R0,2))))

# For environmental noise only
enviro_stats_df_binoms <- enviro_summary_df %>%
  filter(!(name %in% c("max_cases", "duration", "duration_dieout", "peak_time"))) %>%
  rowwise() %>%
  mutate(sum_value = min(sum_value, num_sims)) %>%
  mutate(
    lower_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "lower")$Lower.limit,
    upper_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "upper")$Upper.limit
  ) %>%
  ungroup()

enviro_stats_df_conts <- enviro_summary_df %>%
  filter(name %in% c("max_cases", "duration", "duration_dieout", "peak_time")) %>%
  rowwise() %>%
  mutate(
    max_val = if_else(name == "max_cases", Nh, max_time / 365),
    lower_ci = max(0, mean - 1.96 * sqrt(variance)),
    upper_ci = min(max_val, mean + 1.96 * sqrt(variance))
  ) %>%
  dplyr::select(-max_val) %>%
  ungroup()

enviro_stats_df <- rbind(enviro_stats_df_binoms, enviro_stats_df_conts)
enviro_stats_df$R0_factor = factor(round(enviro_stats_df$R0,2), levels = rev(unique(round(enviro_stats_df$R0,2))))

eps = .Machine$double.eps

# Make comparison dataframes

# Restrict the values of percent difference to assist with visualization
perc_diff_bound = 1

comp_stats_df <- rbind(
  mutate(all_stats_df, type = "all"),
  mutate(enviro_stats_df, type = "enviro")
) %>% 
  mutate(
    .by = c(R0, sigma, name),
    abs_diff = mean[type == "all"] - mean[type == "enviro"],
    perc_diff = if_else(
      mean[type == "enviro"] > 10 * eps, 
      (mean[type == "all"] - mean[type == "enviro"]) / mean[type == "enviro"], # mean_enviro > 0
      if_else(mean[type == "all"] > 10 * eps,  # mean_enviro = 0
              Inf,  # mean_all > 0
              0)    # mean_all < 0
      # (mean[type == "all"] - mean[type == "enviro"]) / mean(mean[type == "enviro"], mean[type == "all"])
    )
  ) %>% 
  mutate(perc_diff = if_else(
    is.infinite(perc_diff), Inf,
    if_else(perc_diff > perc_diff_bound, 
            perc_diff_bound, 
            if_else(perc_diff < -perc_diff_bound,
                    -perc_diff_bound,
                    perc_diff
            ))))

R0s = unique(round(all_stats_df$R0,2))

R0_colors = c(
  # Hotter colors for R0 > 1
  rev(c4a("brewer.yl_or_br", sum(R0s > 1))),
  # Black for R0 == 1
  "black",
  # Cooler colors for R0 < 1
  c4a("scico.oslo", sum(R0s < 1))
)

# [] !!! set up nice plotting to match previous figures

nice_output_labeller = function(output_name) {
  output_label = case_when(
    output_name == "small_outbreak" ~ "Pr(Outbreak > 10 cases)",
    output_name == "big_outbreak" ~ "Pr(Outbreak > 100 cases)",
    output_name == "endemic" ~ "Pr(Transmission lasts 10 years)",
    output_name == "max_cases" ~ "Peak number of cases",
    output_name == "duration" ~ "Outbreak duration",
    output_name == "peak_time" ~ "Peak timing",
    output_name == "duration_dieout" ~ "Outbreak duration (non-endemic case)"
  )
  
  
  return(output_label)
}


# Generic plot function
generic_plot_function <- function(output_name, in_df) {
  in_df %>% 
    ungroup() %>% 
    filter(name == output_name) %>% 
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

# Pr(endemic)
Pr_end_plot <- generic_plot_function("endemic", enviro_stats_df)
ggsave("./figures/no_demo/endemic_prob.png", Pr_end_plot, width = 6.5, height = 3.56525, units = "in")

# Peak cases
Peak_cases_plot <- generic_plot_function("max_cases", enviro_stats_df)
ggsave("./figures/no_demo/peak_cases.png", Peak_cases_plot, width = 6.5, height = 3.56525, units = "in")

# Peak timing
Peak_time_plot <- generic_plot_function("peak_time", enviro_stats_df)
ggsave("./figures/no_demo/peak_time.png", Peak_time_plot, width = 6.5, height = 3.56525, units = "in")

# Small outbreak
Small_outbreak_plot <- generic_plot_function("small_outbreak", enviro_stats_df)
ggsave("./figures/no_demo/small_outbreak_prob.png", Small_outbreak_plot, width = 6.5, height = 3.56525, units = "in")

# Big outbreak
Big_outbreak_plot <- generic_plot_function("big_outbreak", enviro_stats_df)
ggsave("./figures/no_demo/big_outbreak_prob.png", Big_outbreak_plot, width = 6.5, height = 3.56525, units = "in")

# Duration
Duration_plot <- generic_plot_function("duration", enviro_stats_df)
ggsave("./figures/no_demo/duration.png", Duration_plot, width = 6.5, height = 3.56525, units = "in")

# Duration
Duration_dieout_plot <- generic_plot_function("duration_dieout", enviro_stats_df)
ggsave("./figures/duration_dieout.png", Duration_dieout_plot, width = 6.5, height = 3.56525, units = "in")



# Create a "stretched" region at sigma = 0 to show what occurs with no environmental noise
stretch_sigma <- function(in_df) {
  # Step size
  sigma_step = (diff(unique(all_stats_df$sigma)))[1]
  
  # Filter values where sigma = 0
  zero_sigma <- in_df %>% filter(sigma == 0)
  
  # New sigma values
  negative_sigmas = -seq(0.0, 0.25, by = sigma_step)
  
  # Repeat the zero_sigma rows for each negative sigma value
  stretched <- zero_sigma[rep(1:nrow(zero_sigma), each = length(negative_sigmas)), ]
  
  # Assign the new sigma values
  stretched$sigma <- rep(negative_sigmas, times = nrow(zero_sigma))
  
  # Combine original data with stretched region
  bind_rows(in_df, stretched)
}


# Generic plot function
generic_heat_function <- function(output_name, type_name) {
  
  # Set heatmap palette to match unit of output
  palette = case_when(
    output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ "viridis",
    output_name %in% c("duration", "peak_time", "duration_dieout") ~ "plasma",
    output_name %in% c("max_cases") ~ "turbo",
  )
  
  # Label legend appropriately
  leg_label = case_when(
    output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ "Probability",
    output_name %in% c("duration", "peak_time", "duration_dieout") ~ "Time [years]",
    output_name %in% c("max_cases") ~ "Num. cases",
  )
  
  noise_label = case_when(
    type_name == "all" ~  "Demo\n noise only",
    type_name == "enviro" ~ "No\n noise"
  )
  
  
  full_stats_df %>% 
    ungroup() %>% 
    filter(name == output_name, type == type_name) %>% 
    # Stretch out sigma = 0
    stretch_sigma() %>% 
    ggplot(aes(x = sigma, y = R0, z = mean)) +
    geom_raster(aes(fill = mean),
                hjust = unique(diff(full_stats_df$sigma))[2] / 2,
                vjust = unique(diff(full_stats_df$R0))[2] / 2) +
    geom_hline(yintercept = 1, color = "red", lwd = 1) +
    geom_vline(xintercept = 0, color = "black", lwd = 1) +
    # Annotate x-axis for demographic noise
    scale_x_continuous(TeX("Environmental noise strength [$\\sigma$]"),
                       limits = c(-0.25, NA),
                       expand = c(0,0),
                       breaks = c(-0.125, seq(0, 2.0, by = 0.25)),
                       labels = c(noise_label, seq(0, 2.0, by = 0.25))
    ) +
    scale_y_continuous("R0",
                       limits = c(0, NA),
                       expand = c(0,0)) +    
    # color:
    scale_fill_viridis_c(
      name = output_name,
      option = palette
    ) +
    # legend:
    guides(fill = guide_colourbar(
      title = leg_label,
      title.position = "top",
      title.hjust = 0.5,
      barheight = 10,
      show.limits = TRUE
    )) +
    ggtitle(nice_output_labeller(output_name)) +
    theme_half_open(10)
}

full_stats_df <- rbind(
  mutate(all_stats_df, type = "all"),
  mutate(enviro_stats_df, type = "enviro")
)

# Pr(endemic)
Pr_end_heat <- generic_heat_function("endemic", "all")
ggsave("./figures/endemic_prob_heat.png", Pr_end_heat, width = 6.5, height = 3.56525, units = "in")
Pr_end_heat <- generic_heat_function("endemic", "enviro")
ggsave("./figures/no_demo/endemic_prob_heat.png", Pr_end_heat, width = 6.5, height = 3.56525, units = "in")

# Peak cases
Peak_cases_heat <- generic_heat_function("max_cases", "all")
ggsave("./figures/peak_cases_heat.png", Peak_cases_heat, width = 6.5, height = 3.56525, units = "in")
Pr_end_heat <- generic_heat_function("max_cases", "enviro")
ggsave("./figures/no_demo/peak_cases_heat.png", Pr_end_heat, width = 6.5, height = 3.56525, units = "in")

# Peak timing
Peak_time_heat <- generic_heat_function("peak_time", "all")
ggsave("./figures/peak_time_heat.png", Peak_time_heat, width = 6.5, height = 3.56525, units = "in")
Pr_end_heat <- generic_heat_function("peak_time", "enviro")
ggsave("./figures/no_demo/peak_time_heat.png", Pr_end_heat, width = 6.5, height = 3.56525, units = "in")

# Small outbreak
Small_outbreak_heat <- generic_heat_function("small_outbreak", "all")
ggsave("./figures/small_outbreak_prob_heat.png", Small_outbreak_heat, width = 6.5, height = 3.56525, units = "in")
Pr_end_heat <- generic_heat_function("small_outbreak", "enviro")
ggsave("./figures/no_demo/small_outbreak_prob_heat.png", Pr_end_heat, width = 6.5, height = 3.56525, units = "in")

# Big outbreak
Big_outbreak_heat <- generic_heat_function("big_outbreak", "all")
ggsave("./figures/big_outbreak_prob_heat.png", Big_outbreak_heat, width = 6.5, height = 3.56525, units = "in")
Pr_end_heat <- generic_heat_function("big_outbreak", "enviro")
ggsave("./figures/no_demo/big_outbreak_prob_heat.png", Pr_end_heat, width = 6.5, height = 3.56525, units = "in")

# Duration
Duration_heat <- generic_heat_function("duration", "all")
ggsave("./figures/duration_heat.png", Duration_heat, width = 6.5, height = 3.56525, units = "in")
Pr_end_heat <- generic_heat_function("duration", "enviro")
ggsave("./figures/no_demo/duration_heat.png", Pr_end_heat, width = 6.5, height = 3.56525, units = "in")

# Duration dieout
Duration_dieout_heat <- generic_heat_function("duration_dieout", "all")
ggsave("./figures/duration_dieout_heat.png", Duration_dieout_heat, width = 6.5, height = 3.56525, units = "in")
Pr_end_heat <- generic_heat_function("duration_dieout", "enviro")
ggsave("./figures/no_demo/duration_dieout_heat.png", Pr_end_heat, width = 6.5, height = 3.56525, units = "in")


# Comparison plots ----
compare_heat_function <- function(output_name, in_df, type) {
  
  if (type == "abs_diff") {
    leg_label = paste0( "Difference in\n", case_when(
      output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ "Probability",
      output_name %in% c("duration", "peak_time", "duration_dieout") ~ "Time [years]",
      output_name %in% c("max_cases") ~ "Num. cases",
    )
    )
    
    type_label = "absolute difference (w vs. w/o dem. noise)"
    
    # Set heatmap palette to match unit of output
    palette = case_when(
      output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ "viridis",
      output_name %in% c("duration", "peak_time", "duration_dieout") ~ "plasma",
      output_name %in% c("max_cases") ~ "turbo",
    )
    
  } else {
    leg_label = "% difference"
    in_df <- mutate(in_df, perc_diff = 100 * perc_diff)
    palette = "magma"
    
    type_label = "percent difference (w vs. w/o dem. noise)"
  }
  
  
  in_df %>% 
    ungroup() %>% 
    filter(name == output_name,
    ) %>% 
    # Stretch out sigma = 0
    stretch_sigma() %>% # !!! put deterministic stuff here
    ggplot(aes(x = sigma, y = R0, z = !!sym(type))) +
    geom_raster(aes(fill = !!sym(type)),
                hjust = unique(diff(in_df$sigma))[2] / 2,
                vjust = unique(diff(in_df$R0))[2] / 2) +
    # Add a red line for R0 = 1
    geom_hline(yintercept = 1, color = "red", lwd = 1) +
    # Add a black line for "no environmental noise"
    geom_vline(xintercept = 0, color = "black", lwd = 1) +
    # Annotate x-axis for demographic noise
    scale_x_continuous(TeX("Environmental noise strength [$\\sigma$]"),
                       limits = c(-0.25, NA),
                       expand = c(0,0),
                       breaks = c(-0.125, seq(0, 2.0, by = 0.25)),
                       labels = c("No\n noise", seq(0, 2.0, by = 0.25))
    ) +
    scale_y_continuous("R0",
                       limits = c(0, NA),
                       expand = c(0,0)) +    
    # color:
    scale_fill_viridis_c(
      option = palette
    ) +
    # legend:
    guides(fill = guide_colourbar(
      title = leg_label,
      title.position = "top",
      title.hjust = 0.5,
      barheight = 10,
      show.limits = TRUE
    )) +
    ggtitle(paste0(nice_output_labeller(output_name), " ", type_label)) +
    theme_half_open(10)
}

# Pr(endemic)
Pr_end_plot <- compare_heat_function("endemic", comp_stats_df, "abs_diff")
ggsave("./figures/endemic_prob_abs_comparison.png", Pr_end_plot, width = 6.5, height = 3.56525, units = "in")
perc_Pr_end_plot <- compare_heat_function("endemic", comp_stats_df, "perc_diff")
ggsave("./figures/endemic_prob_perc_comparison.png", perc_Pr_end_plot, width = 6.5, height = 3.56525, units = "in")

# Peak cases
Peak_cases_plot <- compare_heat_function("max_cases", comp_stats_df, "abs_diff")
ggsave("./figures/peak_cases_abs_comparison.png", Peak_cases_plot, width = 6.5, height = 3.56525, units = "in")
perc_Peak_cases_plot <- compare_heat_function("max_cases", comp_stats_df, "perc_diff")
ggsave("./figures/peak_cases_perc_comparison.png", perc_Peak_cases_plot, width = 6.5, height = 3.56525, units = "in")

# Peak timing
Peak_time_plot <- compare_heat_function("peak_time", comp_stats_df, "abs_diff")
ggsave("./figures/peak_time_abs_comparison.png", Peak_time_plot, width = 6.5, height = 3.56525, units = "in")
perc_Peak_time_plot <- compare_heat_function("peak_time", comp_stats_df, "perc_diff")
ggsave("./figures/peak_time_perc_comparison.png", perc_Peak_time_plot, width = 6.5, height = 3.56525, units = "in")

# Small outbreak
Small_outbreak_plot <- compare_heat_function("small_outbreak", comp_stats_df, "abs_diff")
ggsave("./figures/small_outbreak_prob_abs_comparison.png", Small_outbreak_plot, width = 6.5, height = 3.56525, units = "in")
perc_Small_outbreak_plot <- compare_heat_function("small_outbreak", comp_stats_df, "perc_diff")
ggsave("./figures/small_outbreak_prob_perc_comparison.png", perc_Small_outbreak_plot, width = 6.5, height = 3.56525, units = "in")

# Big outbreak
Big_outbreak_plot <- compare_heat_function("big_outbreak", comp_stats_df, "abs_diff")
ggsave("./figures/big_outbreak_prob_abs_comparison.png", Big_outbreak_plot, width = 6.5, height = 3.56525, units = "in")
perc_Big_outbreak_plot <- compare_heat_function("big_outbreak", comp_stats_df, "perc_diff")
ggsave("./figures/big_outbreak_prob_perc_comparison.png", perc_Big_outbreak_plot, width = 6.5, height = 3.56525, units = "in")

# Duration
Duration_plot <- compare_heat_function("duration", comp_stats_df, "abs_diff")
ggsave("./figures/duration_abs_comparison.png", Duration_plot, width = 6.5, height = 3.56525, units = "in")
perc_Duration_plot <- compare_heat_function("duration", comp_stats_df, "perc_diff")
ggsave("./figures/duration_perc_comparison.png", perc_Duration_plot, width = 6.5, height = 3.56525, units = "in")

# Duration dieout
Duration_dieout_plot <- compare_heat_function("duration_dieout", comp_stats_df, "abs_diff")
ggsave("./figures/duration_dieout_abs_comparison.png", Duration_dieout_plot, width = 6.5, height = 3.56525, units = "in")
perc_Duration_dieout_plot <- compare_heat_function("duration_dieout", comp_stats_df, "perc_diff")
ggsave("./figures/duration_dieout_perc_comparison.png", perc_Duration_dieout_plot, width = 6.5, height = 3.56525, units = "in")


# All simulations by R0 and sigma
sims_out = read_csv("./data/trajectories_for_grid_plot_no_demo.csv")

All_sims_plot_df <- sims_out %>% # Just use the first 20 simulations
  mutate(R0 = R0_from_Thv_function(Thv)) %>% 
  dplyr::select(-V)

All_sims_plot_df$R0_factor = factor(round(All_sims_plot_df$R0,2), levels = rev((unique(round(All_sims_plot_df$R0,2)))))

endemic_df = All_sims_plot_df %>%   
  group_by(run, R0, sigma) %>% 
  summarise(endemic = H[time == max_time] > 0,
            .groups = "keep")

# Add in whether disease persisted, for coloring
All_sims_plot_df = right_join(All_sims_plot_df, endemic_df)

All_sims_plot <- All_sims_plot_df %>% 
  ggplot(aes(x = time/365, y = H, group = run, color = endemic)) +
  geom_line(alpha = 0.25, lwd = 0.1) +
  facet_grid(rows = vars(R0_factor), cols = vars(round(sigma,2)), scales = 'free') +
  scale_color_manual(values = c("#DF536B", "#61D04F"), breaks = c("TRUE", "FALSE")) +
  scale_x_continuous("Time (years)", limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous("Infected humans", limits = c(0, NA), expand = c(0,0)) +
  guides(color = "none") +
  theme_half_open()

some_sims_plot <- All_sims_plot_df %>% 
  filter(R0_factor %in% c(1.2, 1.1, 1, 0.95, 0.75),
         sigma %in% c(0, 0.2, 0.4, 0.6, 0.8)) %>% 
  ggplot(aes(x = time/365, y = H, group = run, color = endemic)) +
  geom_line(alpha = 0.25, lwd = 0.1) +
  facet_grid(rows = vars(R0_factor), cols = vars(sigma), scales = 'free') +
  scale_color_manual(values = c("#DF536B", "#61D04F"), breaks = c("TRUE", "FALSE")) +
  scale_x_continuous("Time (years)", limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous("Infected humans", limits = c(0, NA), expand = c(0,0)) +
  guides(color = "none") +
  theme_half_open()

ggsave("./figures/no_demo/some_sims.png", some_sims_plot, width = 15, height = 10, units = "in")

ggsave("./figures/no_demo/all_sims.png", All_sims_plot, width = 15, height = 10, units = "in")

# Example trajectory plots ----
# TODO:
# [x] Identify hotspots, coolspots, and nullspots (as R0, sigma pairs)
#    - R0 near 1, sigma between 0.05 and 0.40
#    - (1.05, 0.3), (0.95, 0.3), (2.05, 0.3), 
#    - (0.95, 0.05), (1.05, 0.05), (2.05, 0.05), 
#    - (0.95, 1), (1.05, 1), (2.05, 1)
#    - Make the "0.3" value actually close to the "pit" in the peak cases comparison plot
#       - R0 = 1.375, sigma = 0.55
# [x] Make trajectory comparison plots for chosen R0, sigma values
#    [x] Explaining differences of noise types
#    - 3x3 grid of trajectories with values above
#    - Plot 10-20 sample trajectories (high transparency) and their mean at each time point (low transparency)
#    [x] Explaining "arc" in outbreak duration in the die out cases
#    - Find R0 at max(outbreak_dieout), then the same R0 value with sigma = 1 and 3
#    - Find sigma at max(outbreak_dieout), then the same sigma value with R0 = 1.05 and 3
#    [x] For each of these points, plot distributions of outputs, see how well a Gamma distribution fits.
# [] Set up diverging palettes with set points at the maximum, zero, and minimum
# [] Do a better job with "dieout" cases
# [] - Restrict to cases where duration < 10 years and duration > 0 years
#    - Plot peak cases only when outbreak dies out
# [] In Julia: define outbreak duration as max(time[H > 1]) - min(time[H > 1])
# [] Make more noise difference comparison plots
#    - All noise vs. no environmental (sigma = 0) NB: this is already available
#    - Enviro only vs. deterministic
#    - Dem only vs. deterministic

comparison_trajectories = read.csv("./data/comparison_trajectories.csv") %>% 
  mutate(R0 = R0_from_Thv_function(Thv),
         time = time / 365) %>% 
  group_by(type, R0, sigma, run) %>% 
  mutate(max_H = max(H, na.rm = T),
         max_time = min(time[H == max(H)], na.rm = T)
  ) %>% 
  group_by(type, time, R0, sigma) %>% 
  mutate(mean_H = mean(H, na.rm = T)) %>% 
  ungroup()

comparison_trajectories$R0 = factor(
  comparison_trajectories$R0,
  levels = rev(sort(unique(comparison_trajectories$R0)))
)
comparison_trajectories$type = factor(
  comparison_trajectories$type,
  levels = c("Deterministic", "No_demographic", "All_noise")
)


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
  ggtitle("Points show peak number of cases for each run of the stochastic models") +
  theme_minimal(11) +
  theme_minimal(11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"    
  )


peak_point_plot
ggsave("./figures/peak_cases_comparison.png", peak_point_plot, width = 6.5, height = 4.5, units = "in")

# Peak histogram plot
library(ggh4x)
peak_histogram_plot <- comparison_trajectories %>% 
  filter(type != "Deterministic") %>% 
  dplyr::select(max_time, max_H, type, R0, sigma) %>% 
  unique() %>% 
  ggplot(aes(x = max_H, color = type, fill = type)) +
  # Peak points histogram
  geom_histogram(
    aes(y = after_stat(count)),  # Ensure y-axis represents counts per facet
    # position = "identity",  # Prevent stacking
    position = 'dodge',
    alpha = 0.7,
    bins = 30  # Set fixed bin count (adjustable)
  ) +
  # Means of distributions
  geom_vline(
    data = comparison_trajectories %>% 
      filter(type != "Deterministic") %>% 
      dplyr::select(max_time, max_H, type, R0, sigma) %>% 
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

peak_histogram_plot
ggsave("./figures/peak_cases_comparison_histogram.png", peak_histogram_plot, width = 6.5, height = 4.5, units = "in")

# Explaining "arc" in outbreak duration in the die out cases
#    - Find R0 at max(outbreak_dieout), then the same R0 value with sigma = 1 and 1.5
#    - Find sigma at max(outbreak_dieout), then the same sigma value with R0 = 1.125 and 3
duration_dieout_vals = full_stats_df %>% 
  # ungroup() %>% 
  filter(name == "duration_dieout",
         type == "enviro") %>% 
  pivot_wider(values_from = c(mean)) %>% 
  filter(
    !is.na(duration_dieout),
    duration_dieout > 0
  )

duration_die_out_sigmas = duration_dieout_vals %>% 
  filter(sigma %in% c(1,1.5)) %>% 
  group_by(sigma) %>% 
  filter(duration_dieout == max(duration_dieout, na.rm = T))

duration_die_out_R0s = duration_dieout_vals %>% 
  filter(R0_factor %in% c(1.12,3)) %>% 
  group_by(R0) %>% 
  filter(duration_dieout == max(duration_dieout, na.rm = T))

test_vals = rbind(
  unique(dplyr::select(duration_die_out_sigmas, sigma,R0)),
  unique(dplyr::select(duration_die_out_R0s, sigma, R0))
)

outbreak_duration_df = all_df %>% #enviro_df %>%
  rename(
    max_cases = max_value,
    duration = positive_duration
  ) %>% 
  mutate(
    small_outbreak = ifelse(exceeded_10 == "true", T, F),
    big_outbreak = ifelse(exceeded_100 == "true", T, F),
    endemic = ifelse(positive_at_final == "true", T, F),
    peak_time = max_time / 365, # convert to years
    duration = duration / 365, # convert to years
    duration_dieout = if_else(endemic, NA, 
                              if_else(max_cases <1,
                                      NA,
                                      duration)),
    R0 = R0_from_Thv_function(Thv)
  ) %>% 
  dplyr::select(R0, sigma, run, duration_dieout) %>% 
  filter(
    sigma %in% test_vals$sigma, 
    R0 %in% test_vals$R0,
    !is.na(duration_dieout),
    duration_dieout > 0
  )
outbreak_duration_df$R0_factor = factor(round(outbreak_duration_df$R0,2), levels = rev(unique(round(outbreak_duration_df$R0,2))))

# For each of these points, plot distributions of outputs, see how well a Gamma distribution fits.

# Histogram
duration_dieout_histograms <- outbreak_duration_df %>% 
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

duration_dieout_histograms

# Fit Gamma distribution, perform tests, and generate plot

fitted_data <- outbreak_duration_df %>%
  filter(!is.na(duration_dieout)) %>% 
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
ggplot(fitted_data, aes(x = value)) +
  # Histogram
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.6, color = "black") +
  # Fitted Gamma density curve
  geom_line(data = density_data, aes(x = density_x, y = density_y), color = "red", size = 1) +
  # Gamma distribution fit parameters
  geom_text(data = fitted_data, aes(
    x = Inf, y = Inf,
    label = paste("Shape =", round(shape, 2), "\nRate =", round(rate, 2), "\np-value =", signif(p_value, 3))
  ),
  hjust = 1.1, vjust = 1.1, size = 2, inherit.aes = FALSE
  ) +
  # Means of distributions
  geom_vline(
    data = fitted_data %>% 
      group_by(R0_factor, sigma) %>% 
      summarise(mean_value = mean(value)),
    aes(xintercept = mean_value),
    lwd = 1, lty = 2
  ) +
  # Subplots for each combination of environmental noise strength and R0
  ggh4x::facet_grid2(R0_factor ~ sigma,
             labeller = labeller(.rows = as_labeller(appender_R0,
                                                     default = label_parsed),
                                 .cols = as_labeller(appender_sigma,
                                                     default = label_parsed)),
             scales = "free",
             independent = "all") +
  # scale_x_log10() +
  scale_y_continuous(
    # limits = c(0,1)
  ) +
  labs(
    title = "Histograms with Fitted Gamma Distributions",
    x = "Data Values",
    y = "Density"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12), panel.spacing = unit(1, "lines"))

ggsave("./figures/duration_dieout_fits.png", width = 6.5, height = 4.5, units = "in")