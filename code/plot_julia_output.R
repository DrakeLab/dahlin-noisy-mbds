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
# alphab <- 1 # turn on noise for biting rate
# alphat <- 1 # turn on noise for transmission rate
# alpham <- 1 # turn on noise for mortality rate
# 
# # initial values
# H0 <- 0
# V0 <- 10

# R0 values to test
# R0s <- rev(c(6, 5, 4, 3, 2, 1.25, 1.2, 1.15, 1.1, 1.05, 1, 0.95, 0.75, 0.5, 0))

R0_from_Thv_function <- function(Thv) {
  # Thv <- (R0^2) / ((b^2 * Tvh * Nv) / (Nh * gammah * muv)) 
  sqrt(Thv * b^2 * Tvh * Nv / (Nh * gammah * muv))
}

# Sigma values (environmental noise level in 0-0.3) to test
sigmas <- seq(0, 1, by = 0.05)

# # [] !!! set appropriate resolution
# # Number of batches to run
# num_batches = 100
# # Number of simulations to run
# num_sims = 100
# # Time step size
# deltat = 1
# Final time point
max_time = 3650
# # Time steps 
# time_vec = seq(deltat, max_time, by = deltat)
# # Number of time steps
# num_timesteps = length(time_vec)
# 
# # Dataframes set up ----
# # Set up Variable data frame
# variables_df = expand_grid(
#   tibble(R0 = R0s), 
#   tibble(sigma = sigmas))
# 
# init_pre_df = data.frame(
#   w = as.integer(), 
#   R0 = as.double(), 
#   sigma = as.double(),
#   endemic = as.logical(),
#   small_outbreak = as.logical(),
#   big_outbreak = as.logical(),
#   max_cases = as.double(),
#   duration = as.double()
# )
# 
# init_summary_df = data.frame(
#   R0 = as.double(), 
#   sigma = as.double(),
#   name = as.character(),
#   mean = as.double(),
#   variance = as.double(),
#   sum_value = as.integer()
# )
# 
# 
# # Get data for plotting ----
# 
# # Do you want to run the simulations again (TRUE)?
# # or use old data (FALSE)?
# do_calc = TRUE
# 
# if (do_calc) {
#   # Simulations dataframe
#   init_sims_df = data.frame(
#     w = as.integer(), 
#     batch = as.integer(),
#     R0 = as.double(), 
#     sigma = as.double(),
#     time = as.double(),
#     H = as.double(),
#     V = as.double()
#   )
#   
#   summary_df = init_summary_df
#   
#   # Run simulations ----
#   # Set up parallelization
#   plan(multisession)
#   handlers("cli")
#   
#   for (batch_num in 1:num_batches) {
#     print(paste0("Batch ", batch_num, " / ", num_batches))
#     
#     with_progress({
#       p <- progressor(dim(variables_df)[1])
#       ## Calculate summary statistics across parameter space ----
#       pre_summary_df <- foreach(
#         q = dim(variables_df)[1]:1,
#         .init = init_summary_df,
#         .combine = 'rbind',
#         .options.future = list(seed = TRUE)
#       ) %dofuture% {
#         p()
#         sims_df = init_sims_df
#         # Basic Reproduction Number
#         R0 <- variables_df$R0[q]
#         # Environmental oise level
#         sigma <- variables_df$sigma[q]
#         
#         # calculate Thv based on R0 value
#         Thv <- (R0^2) / ((b^2 * Tvh * Nv) / (Nh * gammah * muv)) 
#         inner_summary_df = init_pre_df
#         for (w in (1:num_sims)) {
#           
#           # set up vectors to store results with value at time=0
#           H_t <- H0
#           V_t <- V0
#           H <- c(H0)
#           V <- c(V0)
#           time = c(0)
#           # calculate random values all at once to avoid repeated calls
#           dW_df = tibble(dW1 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat)),
#                          dW2 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat)),
#                          dW3 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat))
#           )
#           Y = c(0)
#           for (time_index in 1:num_timesteps) {
#             # following equations in paper draft-double check these ****
#             H_s = H_t
#             V_s = V_t
#             
#             dH_det = ((b * Thv / Nh) * V_s * (Nh - H_s) - gammah * H_s)
#             dH_demo = sqrt(((b * Thv / Nh) * V_s * (Nh - H_s) + gammah * H_s))
#             dH_enviro = alphab * b * (Thv / Nh) * V_s * (Nh - H_s)
#             
#             dH = dH_det * deltat + dH_demo * dW_df$dW1[time_index] + 
#               sigma * dH_enviro * dW_df$dW3[time_index]
#             Y_t = H_s + dH
#             Y = append(Y, Y_t)
#             K_t = -min(0, min(Y))
#             H_t = Y_t + K_t 
#             H_t = min(Nh, max(H_s + dH, 0))
#             
#             dV_det = b * (Tvh / Nh) * H_s * (Nv - V_s) - muv * V_s
#             dV_demo = sqrt(b * (Tvh / Nh) * H_s * (Nv - V_s) + muv * V_s)
#             dV_enviro = alphab * b * (Tvh / Nh) * H_s * (Nv - V_s) + alphat * b * (Tvh / Nh) * H_s * (Nv - V_s) + alpham * muv * V_s
#             
#             dV = dV_det * deltat + dV_demo * dW_df$dW2[time_index] + 
#               sigma * dV_enviro * dW_df$dW3[time_index]
#             V_t = min(Nv, max(V_s + dV, 0))
#             
#             # Add the values to vector
#             H <- append(H, H_t)
#             V <- append(V, V_t)
#             time <- append(time, time_vec[time_index])
#             
#             # if both H_t and V_t are 0, end the simulation.
#             if (H_t == 0 & V_t == 0) {
#               H <- append(H, 0)
#               V <- append(V, 0)
#               time <- append(time, max_time)
#               break
#             }
#           }
#           
#           # Calculate summary statistics
#           pre_summary_row = data.frame(
#             w = w, 
#             batch = batch_num,
#             R0 = R0, 
#             sigma = sigma,
#             time = time,
#             H = H,
#             V = V
#           ) %>% 
#             # Calculate statistics for each simulation
#             summarise(
#               .by = c(R0, sigma, w, batch),
#               endemic = H[time == max_time] > 0, #
#               small_outbreak = any(H > 10),
#               big_outbreak = any(H > 100),
#               max_cases = max(H),
#               duration = ifelse(
#                 max_cases > 1, # If an outbreak actually occurred...
#                 max(time[H > 1 | V > 1]) - min(time[H > 1 | V > 1]), # Measure how long cases were above 1
#                 0 # Otherwise the duration is zero
#               )
#             )
#           
#           inner_summary_df = rbind(inner_summary_df, pre_summary_row)
#         }
#         inner_summary_df
#       }
#     })
#     pre_summary_df = rbind(pre_summary_df, pre_summary_df)
#   }
#   plan(sequential)
#   gc()
#   
#   summary_df = pre_summary_df %>% 
#     # Summarise for each combination of R0 and sigma
#     pivot_longer(cols = endemic:duration) %>% 
#     mutate(
#       .by = c(R0, sigma, name),
#       mean = mean(value),
#       variance = var(value),
#       sum_value = as.integer(sum(value)) # needed for Clopper-Pearon confidence intervals
#     ) %>% 
#     select(-c(w, value)) %>% distinct() 
#   # summary_rows
#   
#   summary_df$R0_factor = factor(summary_df$R0, levels = rev(R0s))
#   
#   summary_df_binoms <- summary_df %>%
#     filter(!(name %in% c("max_cases", "duration"))) %>% 
#     rowwise() %>% 
#     mutate(sum_value = min(sum_value, num_sims)) %>% 
#     mutate(
#       lower_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "lower")$Lower.limit,
#       upper_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "upper")$Upper.limit
#     ) %>% 
#     ungroup()
#   
#   summary_df_conts <- summary_df %>%
#     filter(name %in% c("max_cases", "duration")) %>% 
#     rowwise() %>% 
#     mutate(
#       max_val = if_else(name == "max_cases", Nh, max_time),
#       lower_ci = max(0, mean - 1.96 * sqrt(variance)),
#       upper_ci = min(max_val, mean + 1.96 * sqrt(variance))
#     ) %>% 
#     select(-max_val) %>% 
#     ungroup()
#   
#   final_summary_df <- rbind(summary_df_binoms, summary_df_conts)
#   write_rds(final_summary_df, "./data/summary_data.rds")
#   
#   ## Calculate some individual simulations to plot trajectories ----
#   plan(multisession)
#   handlers("cli")
#   with_progress({
#     p <- progressor(dim(variables_df)[1])
#     
#     sims_out <- foreach(
#       q = 1:dim(variables_df)[1],
#       .init = init_sims_df,
#       .combine = 'rbind',
#       .options.future = list(seed = TRUE)
#     ) %dofuture% {
#       p()
#       sims_df = init_sims_df
#       # Basic Reproduction Number
#       R0 <- variables_df$R0[q]
#       # Environmental oise level
#       sigma <- variables_df$sigma[q]
#       
#       # calculate Thv based on R0 value
#       Thv <- (R0^2) / ((b^2 * Tvh * Nv) / (Nh * gammah * muv)) 
#       pre_summary_df = init_pre_df
#       for (w in (1:50)) {
#         # set up vectors to store results with value at time=0
#         H_t <- H0
#         V_t <- V0
#         H <- c(H0)
#         V <- c(V0)
#         time = c(0)
#         # calculate random values all at once to avoid repeated calls
#         dW_df = tibble(dW1 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat)),
#                        dW2 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat)),
#                        dW3 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat))
#         )
#         
#         for (time_index in 1:length(time_vec)) {
#           # following equations in paper draft-double check these ****
#           H_s = H_t
#           V_s = V_t
#           
#           dH_det = ((b * Thv / Nh) * V_s * (Nh - H_s) - gammah * H_s)
#           dH_demo = sqrt(((b * Thv / Nh) * V_s * (Nh - H_s) + gammah * H_s))
#           dH_enviro = alphab * b * (Thv / Nh) * V_s * (Nh - H_s)
#           
#           dH = dH_det * deltat + dH_demo * dW_df$dW1[time_index] + 
#             sigma * dH_enviro * dW_df$dW3[time_index]
#           H_t = min(Nh, max(H_s + dH, 0))
#           
#           dV_det = b * (Tvh / Nh) * H_s * (Nv - V_s) - muv * V_s
#           dV_demo = sqrt(b * (Tvh / Nh) * H_s * (Nv - V_s) + muv * V_s)
#           dV_enviro = alphab * b * (Tvh / Nh) * H_s * (Nv - V_s) + alphat * b * (Tvh / Nh) * H_s * (Nv - V_s) + alpham * muv * V_s
#           
#           dV = dV_det * deltat + dV_demo * dW_df$dW2[time_index] + 
#             sigma * dV_enviro * dW_df$dW3[time_index]
#           V_t = min(Nv, max(V_s + dV, 0))
#           
#           # Add the values to vector
#           H <- append(H, H_t)
#           V <- append(V, V_t)
#           time <- append(time, time_vec[time_index])
#           
#           # if both H_t and V_t are 0, end the simulation.
#           if (H_t == 0 & V_t == 0) {
#             H <- append(H, 0)
#             V <- append(V, 0)
#             time <- append(time, max_time)
#             break
#           }
#         }
#         
#         # Calculate summary statistics
#         sim_row = data.frame(
#           w = w, 
#           R0 = R0, 
#           sigma = sigma,
#           time = time,
#           H = H,
#           V = V
#         )
#         sims_df = rbind(sims_df, sim_row)
#       }
#       sims_df
#     }
#   })
#   
#   plan(sequential)
#   gc()
#   sims_out$R0_factor = factor(sims_out$R0, levels = rev(R0s))
#   write_rds(sims_out %>% filter(w < 101), "./data/some_sims.rds")
#   
#   
# } else {
#   summary_df = read_rds("./data/summary_data.rds")
#   sims_out = read_rds("./data/some_sims.rds")
# }

# Plots ----


read_df = read.csv("./data/collect_all_outputs.csv") 
summary_df = read_df %>%
  rename(
    max_cases = max_value,
    duration = positive_duration
  ) %>% 
  mutate(
    small_outbreak = ifelse(exceeded_10 == "true", T, F),
    big_outbreak = ifelse(exceeded_100 == "true", T, F),
    endemic = ifelse(positive_at_final == "true", T, F),
    duration = duration / 365, # convert to years
    R0 = R0_from_Thv_function(Thv),
    .keep = "unused"
  ) %>% 
  # Summarise for each combination of R0 and sigma
  pivot_longer(cols = max_cases:endemic) %>%
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
  filter(!(name %in% c("max_cases", "duration"))) %>%
  rowwise() %>%
  mutate(sum_value = min(sum_value, num_sims)) %>%
  mutate(
    lower_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "lower")$Lower.limit,
    upper_ci = clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "upper")$Upper.limit
  ) %>%
  ungroup()

stats_df_conts <- summary_df %>%
  filter(name %in% c("max_cases", "duration")) %>%
  rowwise() %>%
  mutate(
    max_val = if_else(name == "max_cases", Nv, max_time),
    lower_ci = max(0, mean - 1.96 * sqrt(variance)),
    upper_ci = min(max_val, mean + 1.96 * sqrt(variance))
  ) %>%
  select(-max_val) %>%
  ungroup()

stats_df <- rbind(stats_df_binoms, stats_df_conts)
stats_df$R0_factor = factor(round(stats_df$R0,2), levels = rev(unique(round(stats_df$R0,2))))

R0_colors = c(
  # Hotter colors for R0 > 1
  rev(c4a("brewer.yl_or_br", 16)),
  # Black for R0 == 1
  "black",
  # Cooler colors for R0 < 1
  c4a("scico.oslo", 4)
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
    geom_contour_filled(bins = 20) +
    geom_hline(yintercept = 1, color = "red", lwd = 1) +
    # Confidence intervals
    scale_x_continuous(TeX("Environmental noise strength [$\\sigma$]"),
                       limits = c(0,NA),
                       expand = c(0,0)) +
    scale_y_continuous("R0",
                       limits = c(0,NA),
                       expand = c(0,0)) +    
    # color:
    scale_fill_viridis_d(
      name = output_name,
      option = "plasma"
    ) +
    # legend:
    guides(fill = guide_coloursteps(
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
sims_out = read_csv("./code/julia//trajectories_for_grid_plot.csv")

All_sims_plot_df <- sims_out %>% # Just use the first 20 simulations
  mutate(R0 = R0_from_Thv_function(Thv)) %>% 
  select(-V)

All_sims_plot_df$R0_factor = factor(round(All_sims_plot_df$R0,2), levels = (unique(round(All_sims_plot_df$R0,2))))

endemic_df = All_sims_plot_df %>%   
  group_by(run, R0, sigma) %>% 
  summarise(endemic = H[time == max_time] > 0,
            .groups = "keep")

# Add in whether disease persisted, for coloring
All_sims_plot_df = right_join(All_sims_plot_df, endemic_df)

All_sims_plot <- All_sims_plot_df %>% 
  ggplot(aes(x = time/365, y = H, group = run, color = endemic)) +
  geom_line(alpha = 0.25, lwd = 0.1) +
  facet_grid(rows = vars(R0_factor), cols = vars(sigma), scales = 'free') +
  scale_color_manual(values = c("#DF536B", "#61D04F"), breaks = c("TRUE", "FALSE")) +
  scale_x_continuous("Time (years)", limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous("Infected humans", limits = c(0, NA), expand = c(0,0)) +
  guides(color = "none") +
  theme_half_open(4)

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
  theme_half_open(4)

ggsave("./figures/some_sims.png", some_sims_plot, width = 6.5, height = 3.56525, units = "in")

ggsave("./figures/all_sims.png", All_sims_plot, width = 6.5, height = 3.56525, units = "in")