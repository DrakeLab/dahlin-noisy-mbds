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

# Labeling functions
appender_R0 <- function(string) TeX(paste("$R_0 = $", string))
appender_sigma <- function(string) TeX(paste("$\\sigma = $", string))

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

# Load data from Julia output ----
# Deterministic
det_df = read_csv("./data/collect_all_outputs_det.csv.gz") 
# Environmental noise only
enviro_df = read_csv("./data/collect_all_outputs_no_demo.csv.gz") 
# Demographic + Environmental noise
all_df = read_csv("./data/collect_all_outputs.csv.gz")

# # Set up parallel processing
# n_cores <- parallel::detectCores() - 1
# cluster <- new_cluster(n_cores)
# cluster_library(cluster, "dplyr")
# # Get all functions in the global environment
# # all_stuff <- ls(envir = globalenv(), pattern = "^[a-zA-Z]")  
# cluster_copy(cluster, c("R0_from_Thv_function", "b", "Tvh", "Nv", "Nh", "gammah", "muv"))

# Summarize key statistics ----
# calculate mean, variance, and number of positives
all_df_modified <- all_df %>%
  # partition(cluster) %>% 
  rename(
    max_cases = max_value,
    duration = positive_duration,
    small_outbreak = exceeded_10,
    big_outbreak = exceeded_100,
    endemic = positive_at_final
  ) %>% 
  mutate(
    peak_time = max_time / 365, # convert to years
    duration = duration / 365, # convert to years
    duration_10 = if_else(small_outbreak & !endemic, duration, NA),
    duration_100 = if_else(big_outbreak & !endemic, duration, NA),
    duration_dieout = if_else(endemic, NA, 
                              if_else(max_cases <1,
                                      NA,
                                      duration)),
    R0 = R0_from_Thv_function(Thv)
  ) %>% 
  # collect() %>% 
  dplyr::select(-c(max_time, Thv)) %>% 
  # Summarise for each combination of R0 and sigma
  pivot_longer(cols = max_cases:duration_dieout)
rm(all_df); gc()
write_rds(all_df_modified, "./data/all_modified.rds")

# start_time <- Sys.time()
all_summary_df = all_df_modified %>%
  # partition(cluster) %>% 
  mutate(
    .by = c(R0, sigma, name),
    mean = mean(value, na.rm = TRUE),
    variance = var(value, na.rm = TRUE),
    sum_value = as.integer(sum(value, na.rm = TRUE)) # needed for Clopper-Pearon confidence intervals
  ) %>%
  # collect() %>% 
  dplyr::select(-c(run, value)) %>% distinct()
# end_time <- Sys.time()
# print(end_time - start_time)

enviro_df_modified <- enviro_df %>%
  rename(
    max_cases = max_value,
    duration = positive_duration,
    small_outbreak = exceeded_10,
    big_outbreak = exceeded_100,
    endemic = positive_at_final
  ) %>% 
  mutate(
    peak_time = max_time / 365, # convert to years
    duration = duration / 365, # convert to years
    duration_10 = if_else(small_outbreak & !endemic, duration, NA),
    duration_100 = if_else(big_outbreak & !endemic, duration, NA),
    duration_dieout = if_else(endemic, NA, 
                              if_else(max_cases <1,
                                      NA,
                                      duration)),
    R0 = R0_from_Thv_function(Thv)
  ) %>% 
  dplyr::select(-c(max_time, Thv)) %>% 
  # Summarise for each combination of R0 and sigma
  pivot_longer(cols = max_cases:duration_dieout) 
rm(enviro_df); gc()
write_rds(enviro_df_modified, "./data/enviro_modified.rds")

enviro_summary_df = enviro_df_modified %>%
  mutate(
    .by = c(R0, sigma, name),
    mean = mean(value, na.rm = TRUE),
    variance = var(value, na.rm = TRUE),
    sum_value = as.integer(sum(value, na.rm = TRUE)) # needed for Clopper-Pearon confidence intervals
  ) %>%
  dplyr::select(-c(run, value)) %>% distinct()

# Process statistics ----
# Calculate confidence intervals
num_sims = max(enviro_df_modified$run)
# For both noise types
all_stats_df_binoms <- all_summary_df %>%
  filter(!(name %in% c("max_cases", "duration", "duration_dieout", "peak_time", "duration_10", "duration_100"))) %>%
  rowwise() %>%
  mutate(sum_value = min(sum_value, num_sims)) %>%
  mutate(
    lower_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "lower")$Lower.limit,
    upper_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "upper")$Upper.limit
  ) %>%
  ungroup()

all_stats_df_conts <- all_summary_df %>%
  filter(name %in% c("max_cases", "duration", "duration_dieout", "peak_time", "duration_10", "duration_100")) %>%
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
write_rds(all_stats_df, "./data/all_stats.rds")

# For environmental noise only
enviro_stats_df_binoms <- enviro_summary_df %>%
  filter(!(name %in% c("max_cases", "duration", "duration_dieout", "peak_time", "duration_10", "duration_100"))) %>%
  rowwise() %>%
  mutate(sum_value = min(sum_value, num_sims)) %>%
  mutate(
    lower_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "lower")$Lower.limit,
    upper_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "upper")$Upper.limit
  ) %>%
  ungroup()

enviro_stats_df_conts <- enviro_summary_df %>%
  filter(name %in% c("max_cases", "duration", "duration_dieout", "peak_time", "duration_10", "duration_100")) %>%
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

write_rds(enviro_stats_df, "./data/enviro_stats.rds")

# Make comparison dataframes ----

# Restrict the values of percent difference to assist with visualization
perc_diff_bound = 1
eps = .Machine$double.eps

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

write_rds(comp_stats_df, "./data/comp_stats.rds")

# Trajectories comparisons ----
comparison_trajectories = read_csv("./data/comparison_trajectories.csv") %>% 
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

# Save
write_rds(comparison_trajectories, "./data/comp_trajectories.rds")

# All simulations ----
sims_out = read_csv("./data/trajectories_for_grid_plot_no_demo.csv.gz")

All_sims_plot_df <- sims_out %>% # Just use the first 20 simulations
  mutate(R0 = R0_from_Thv_function(Thv)) %>% 
  dplyr::select(-V)

All_sims_plot_df$R0_factor = factor(round(All_sims_plot_df$R0,2), levels = rev((unique(round(All_sims_plot_df$R0,2)))))

write_rds(All_sims_plot_df, "./data/all_sims.rds")