library(tidyverse)
library(GenBinomApps) # to get Clopper-Pearson confidence interval function
library(cols4all)
library(future)
library(doFuture)
library(foreach)
library(progressr) # progress bars for long parallel computations
library(cowplot)
library(latex2exp)

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
read_df = read.csv("./data/collect_all_outputs.csv") 

# Summarize key statistics: calculate mean, variance, and number of positives
summary_df = read_df %>%
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
    R0 = R0_from_Thv_function(Thv),
    .keep = "unused"
  ) %>% 
  # Summarise for each combination of R0 and sigma
  pivot_longer(cols = max_cases:peak_time) %>%
  mutate(
    .by = c(R0, sigma, name),
    mean = mean(value),
    variance = var(value),
    sum_value = as.integer(sum(value)) # needed for Clopper-Pearon confidence intervals
  ) %>%
  select(-c(run, value)) %>% distinct()

num_sims = max(read_df$run)
# Process statistics
stats_df_binoms <- summary_df %>%
  filter(!(name %in% c("max_cases", "duration", "peak_time"))) %>%
  rowwise() %>%
  mutate(sum_value = min(sum_value, num_sims)) %>%
  mutate(
    lower_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "lower")$Lower.limit,
    upper_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "upper")$Upper.limit
  ) %>%
  ungroup()

stats_df_conts <- summary_df %>%
  filter(name %in% c("max_cases", "duration", "peak_time")) %>%
  rowwise() %>%
  mutate(
    max_val = if_else(name == "max_cases", Nh, max_time / 365),
    lower_ci = max(0, mean - 1.96 * sqrt(variance)),
    upper_ci = min(max_val, mean + 1.96 * sqrt(variance))
  ) %>%
  select(-max_val) %>%
  ungroup()

stats_df <- rbind(stats_df_binoms, stats_df_conts)
stats_df$R0_factor = factor(round(stats_df$R0,2), levels = rev(unique(round(stats_df$R0,2))))

R0s = unique(round(stats_df$R0,2))

R0_colors = c(
  # Hotter colors for R0 > 1
  rev(c4a("brewer.yl_or_br", sum(R0s > 1))),
  # Black for R0 == 1
  "black",
  # Cooler colors for R0 < 1
  c4a("scico.oslo", sum(R0s < 1))
)

# [] !!! set up nice plotting to match previous figures

# Generic plot function
generic_plot_function <- function(output_name) {
  stats_df %>% 
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
Pr_end_plot <- generic_plot_function("endemic")
ggsave("./figures/endemic_prob.png", Pr_end_plot, width = 6.5, height = 3.56525, units = "in")

# Peak cases
Peak_cases_plot <- generic_plot_function("max_cases")
ggsave("./figures/peak_cases.png", Peak_cases_plot, width = 6.5, height = 3.56525, units = "in")

# Peak timing
Peak_time_plot <- generic_plot_function("peak_time")
ggsave("./figures/peak_time.png", Peak_time_plot, width = 6.5, height = 3.56525, units = "in")

# Small outbreak
Small_outbreak_plot <- generic_plot_function("small_outbreak")
ggsave("./figures/small_outbreak_prob.png", Small_outbreak_plot, width = 6.5, height = 3.56525, units = "in")

# Big outbreak
Big_outbreak_plot <- generic_plot_function("big_outbreak")
ggsave("./figures/big_outbreak_prob.png", Big_outbreak_plot, width = 6.5, height = 3.56525, units = "in")

# Duration
Duration_plot <- generic_plot_function("duration")
ggsave("./figures/duration.png", Duration_plot, width = 6.5, height = 3.56525, units = "in")

# Generic plot function
generic_heat_function <- function(output_name) {
  stats_df %>% 
    ungroup() %>% 
    filter(name == output_name) %>% 
    ggplot(aes(x = sigma, y = R0, z = mean)) +
    geom_tile(aes(fill = mean)) +
    # geom_contour_filled(bins = 20) +
    geom_hline(yintercept = 1, color = "red", lwd = 1) +
    # Confidence intervals
    scale_x_continuous(TeX("Environmental noise strength [$\\sigma$]"),
                       limits = c(0,NA),
                       expand = c(0,0)) +
    scale_y_continuous("R0",
                       limits = c(0,NA),
                       expand = c(0,0)) +    
    # color:
    scale_fill_viridis_c(
      name = output_name,
      option = "plasma"
    ) +
    # legend:
    guides(fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = 24,
      show.limits = TRUE
    )) +
    ggtitle(output_name) +
    theme_half_open()
}

# Pr(endemic)
Pr_end_heat <- generic_heat_function("endemic")
ggsave("./figures/endemic_prob_heat.png", Pr_end_heat, width = 6.5, height = 3.56525, units = "in")

# Peak cases
Peak_cases_heat <- generic_heat_function("max_cases")
ggsave("./figures/peak_cases_heat.png", Peak_cases_heat, width = 6.5, height = 3.56525, units = "in")

# Peak timing
Peak_time_heat <- generic_heat_function("peak_time")
ggsave("./figures/peak_time_heat.png", Peak_time_heat, width = 6.5, height = 3.56525, units = "in")

# Small outbreak
Small_outbreak_heat <- generic_heat_function("small_outbreak")
ggsave("./figures/small_outbreak_prob_heat.png", Small_outbreak_heat, width = 6.5, height = 3.56525, units = "in")

# Big outbreak
Big_outbreak_heat <- generic_heat_function("big_outbreak")
ggsave("./figures/big_outbreak_prob_heat.png", Big_outbreak_heat, width = 6.5, height = 3.56525, units = "in")

# Duration
Duration_heat <- generic_heat_function("duration")
ggsave("./figures/duration_heat.png", Duration_heat, width = 6.5, height = 3.56525, units = "in")














# All simulations by R0 and sigma
sims_out = read_csv("./data/trajectories_for_grid_plot.csv")

All_sims_plot_df <- sims_out %>% # Just use the first 20 simulations
  mutate(R0 = R0_from_Thv_function(Thv)) %>% 
  select(-V)

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

ggsave("./figures/some_sims.png", some_sims_plot, width = 15, height = 10, units = "in")

ggsave("./figures/all_sims.png", All_sims_plot, width = 15, height = 10, units = "in")