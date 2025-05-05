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

# Summarize key statistics ----

all_summary_df <- all_df %>%
  # rename quantities to match plots
  mutate(name = case_when(
    name == "max_value" ~ "max_cases",
    # name == "positive_duration" ~ "duration",
    name == "exceeded_10" ~ "small_outbreak",
    name == "exceeded_100" ~ "big_outbreak",
    name == "positive_at_final" ~ "endemic",
    TRUE ~ name
  )) %>% 
  mutate(R0 = round(R0_from_Thv_function(Thv), 3)) %>% 
  mutate(sigma = round(sigma, 3)) %>% 
  # collect() %>% 
  dplyr::select(-Thv)

enviro_summary_df <- enviro_df %>%
  # rename quantities to match plots
  mutate(name = case_when(
    name == "max_value" ~ "max_cases",
    # name == "positive_duration" ~ "duration",
    name == "exceeded_10" ~ "small_outbreak",
    name == "exceeded_100" ~ "big_outbreak",
    name == "positive_at_final" ~ "endemic",
    TRUE ~ name
  )) %>% 
  mutate(R0 = round(R0_from_Thv_function(Thv), 3)) %>% 
  mutate(sigma = round(sigma, 3)) %>% 
  # collect() %>% 
  dplyr::select(-Thv)


# Process statistics ----
# Calculate confidence intervals
# For both noise types

all_stats_df <- all_summary_df %>% 
  pivot_wider(names_from = statistic) %>% 
  rowwise() %>%
  mutate(
    max_val = if_else(name == "max_cases", Nh, max_time / 365),
    lower_ci = max(0, mean - 0.674 * sqrt(variance)),
    upper_ci = if_else(name %in% c("small_outbreak", "big_outbreak", "endemic", "zero_cases"), 
                       min(1, mean + 0.674 * sqrt(variance)),
                       min(max_val, mean + 0.674 * sqrt(variance))
    )
  ) %>%
  dplyr::select(-max_val) %>%
  ungroup()

all_stats_df$R0_factor = factor(round(all_stats_df$R0,3), levels = rev(unique(round(all_stats_df$R0,3))))
write_rds(all_stats_df, "./data/all_stats.rds")

# For environmental noise only
enviro_stats_df <- enviro_summary_df %>% 
  pivot_wider(names_from = statistic) %>% 
  rowwise() %>%
  mutate(
    max_val = if_else(name == "max_cases", Nh, max_time / 365),
    lower_ci = max(0, mean - 0.674 * sqrt(variance)),
    upper_ci = if_else(name %in% c("small_outbreak", "big_outbreak", "endemic", "zero_cases"), 
                       min(1, mean + 0.674 * sqrt(variance)),
                       min(max_val, mean + 0.674 * sqrt(variance))
    )
  ) %>%
  dplyr::select(-max_val) %>%
  ungroup()

enviro_stats_df$R0_factor = factor(round(enviro_stats_df$R0,3), levels = rev(unique(round(enviro_stats_df$R0,3))))

write_rds(enviro_stats_df, "./data/enviro_stats.rds")

# Make comparison dataframes ----

# Restrict the values of percent difference to assist with visualization
perc_diff_bound = 1
eps = .Machine$double.eps

comp_stats_df <- rbind(
  mutate(all_stats_df, type = "all"),
  mutate(enviro_stats_df, type = "enviro")
) %>% 
  ungroup() %>% 
  mutate(
    .by = c(R0, sigma, name),
    abs_diff = ifelse((is.na(mean[type == "all"]) || is.na(mean[type == "enviro"])),
                      NA,
                      mean[type == "all"] - mean[type == "enviro"]),
    perc_diff = ifelse((is.na(mean[type == "all"]) || is.na(mean[type == "enviro"])),
                        NA,
                        if_else(
      mean[type == "enviro"] > 10 * eps, 
      (mean[type == "all"] - mean[type == "enviro"]) / mean[type == "enviro"], # mean_enviro > 0
      if_else(mean[type == "all"] > 10 * eps,  # mean_enviro = 0
              Inf,  # mean_all > 0
              0)    # mean_all < 0
      # (mean[type == "all"] - mean[type == "enviro"]) / mean(mean[type == "enviro"], mean[type == "all"])
    ))
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
sims_out = read_csv("./data/trajectories_for_grid_plot.csv.gz")

All_sims_plot_df <- sims_out %>% 
  mutate(R0 = R0_from_Thv_function(Thv)) %>% 
  filter(
    H>=0, V>=0,
    H<=Nh, V<=Nv
  ) %>%
  # dplyr::select(-V) %>% 
  group_by(run, R0, sigma) %>% 
  distinct() %>% 
  mutate(endemic = (max(time) == max_time && H[time == max_time] > 1)) %>% 
  ungroup()



All_sims_plot_df$R0_factor = factor(round(All_sims_plot_df$R0,3), levels = rev((unique(round(All_sims_plot_df$R0,3)))))
All_sims_plot_df$sigma_factor = factor(round(All_sims_plot_df$sigma,3), levels = unique(round(All_sims_plot_df$sigma,3)))

write_rds(All_sims_plot_df, "./data/all_sims.rds")

# Duration vs. intensity direct comparisons ----
enviro_outbreak_duration_df = read_csv("./data/dur_peak_no_demo.csv.gz")  %>% 
  mutate(R0 = R0_from_Thv_function(Thv)) %>% 
  filter(
    # round(sigma, 3) %in% c(1, 1.5), #test_vals$sigma,
    # round(R0, 3) %in% c(1.375, 3)#test_vals$R0
  ) %>% 
  dplyr::select(sigma, R0, max_value, max_time, duration_dieout) %>% 
  mutate(type = "Environmental noise only")
enviro_outbreak_duration_df$R0_factor = factor(
  round(enviro_outbreak_duration_df$R0, 3), 
  levels = (unique(round(enviro_outbreak_duration_df$R0, 3))))

all_outbreak_duration_df = read_csv("./data/dur_peak.csv.gz")  %>% 
  mutate(R0 = R0_from_Thv_function(Thv)) %>%
  filter(
    # round(sigma, 3) %in% c(1, 1.5), #test_vals$sigma,
    # round(R0, 3) %in% c(1.375, 3)#test_vals$R0
  ) %>% 
  dplyr::select(sigma, R0, max_value, max_time, duration_dieout) %>% 
  mutate(type = "All noise")
all_outbreak_duration_df$R0_factor = factor(
  round(all_outbreak_duration_df$R0, 3), 
  levels = (unique(round(all_outbreak_duration_df$R0, 3))))

outbreak_duration_df = bind_rows(
  all_outbreak_duration_df,
  enviro_outbreak_duration_df
)

write_rds(outbreak_duration_df, "./data/peak_v_duration_sims.rds")

