library(tidyverse)
library(GenBinomApps) # to get Clopper-Pearson confidence interval function
library(cols4all)
library(future)
library(foreach)
library(progressr) # progress bars for long parallel computations

# Parameters ---- 
b <- 0.3 # biting rate
Tvh <- 0.5 # vector transmission-probability susceptible vector being infected after contact with infectious host
Nh <- 10000 # number of hosts in system
Nv <- 100000 # number of vectors in system
muv <- 0.1 # vector mortality rate
gammah <- 0.1 # host recovery rate
alphab <- 1 # turn on noise for biting rate
alphat <- 1 # turn on noise for transmission rate
alpham <- 1 # turn on noise for mortality rate

# initial values
H0 <- 0
V0 <- 10

# R0 values to test
R0s = c(6, 5, 4, 3, 2, 1.25, 1.2, 1.15, 1.1, 1.05, 1, 0.95, 0.75, 0.5, 0)

# Sigma values (environmental noise level in 0-0.3) to test
sigmas <- seq(0, 0.3, by = 0.05)

# [] !!! set appropriate resolution
# Number of simulations to run
num_sims = 10
# Time step size
deltat = .1
# Final time point
max_time = 3650
# Time steps 
time_vec = seq(deltat, max_time, by = deltat)
# Number of time steps
num_timesteps = length(time_vec)


# Dataframes set up ----
# Set up Variable data frame
variables_df = expand_grid(tibble(R0 = R0s), tibble(sigma = sigmas))

# Simulations dataframe
init_sims_df = data.frame(
  w = as.integer(), 
  R0 = as.double(), 
  sigma = as.double(),
  time = as.double(),
  H = as.double(),
  V = as.double()
)

# Set up end dataframe
init_dataframe <- data.frame(
  R0 = as.double(), 
  sigma = as.double(), 
  # probability of endemic disease-lasting until end of simulation
  endemic = as.logical(), 
  # peak number of cases
  max_cases = as.double(),
  # probability of outbreaks greater than 10 hosts
  small_outbreak = as.logical(), 
  # probability of outbreaks greater than 100 hosts
  big_outbreak = as.logical(),
  # outbreak duration
  duration = as.double()
)
# Final dataframe
results = init_dataframe

# Set up final summary dataframe
# [] Parallelize the loop
# [] !!! set up progress bar

plan(multisession)
handlers(global = TRUE)
sims_out <- foreach(
  q = 1:dim(variables_df)[1],
  .init = init_sims_df,
  .combine = 'rbind',
  .options.future = list(seed = TRUE)
) %dofuture% {
  sims_df = init_sims_df
  # Basic Reproduction Number
  R0 <- variables_df$R0[q]
  # Environmental oise level
  sigma <- variables_df$sigma[q]
  
  # calculate Thv based on R0 value
  Thv <- (R0^2) / ((b^2 * Tvh * Nv) / (Nh * gammah * muv)) 
  
  
  sim_summary =  cbind(data.frame(w = as.integer()), init_dataframe)
  
  # set up vectors to store results with value at time=0
  H <- c(H0)
  V <- c(V0)
  time <- c(0)
  
  H_t <- H0
  V_t <- V0
  
  for (w in (1:num_sims)) {
    H_t <- H0
    V_t <- V0
    H <- c(H0)
    V <- c(V0)
    time = c(0)
    # calculate random values all at once to avoid repeated calls
    dW_df = tibble(dW1 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat)),
                   dW2 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat)),
                   dW3 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat))
    )
    
    for (time_index in 1:length(time_vec)) {
      # following equations in paper draft-double check these ****
      dW1 = dW_df$dW1[time_index]
      dW2 = dW_df$dW2[time_index]
      dW3 = dW_df$dW3[time_index]
      
      H_s <- H_t
      V_s <- V_t
      
      dH_det = ((b * Thv / Nh) * V_s * (Nh - H_s) - gammah * H_s)
      dH_demo = sqrt(((b * Thv / Nh) * V_s * (Nh - H_s) + gammah * H_s))
      dH_enviro = alphab * (Thv / Nh) * V_s * (Nh - H_s)
      
      dH <- dH_det * deltat + dH_demo * dW1 + sigma * dH_enviro * dW3
      H_t = max(H_s + dH, 0)
      
      dV_det = ((b * Thv / Nh) * V_s * (Nh - H_s) - gammah * H_s)
      dV_demo = sqrt(((b * Thv / Nh) * V_s * (Nh - H_s) + gammah * H_s))
      dV_enviro = alphab * (Thv / Nh) * V_s * (Nh - H_s)
      
      dV <- dV_det * deltat + dV_demo * dW2 + sigma * dV_enviro * dW3
      V_t = max(V_s + dV, 0)
      
      # Add the values to vector
      H <- append(H, H_t)
      V <- append(V, V_t)
      time <- append(time, time_vec[time_index])
      
      # if both H_t and V_t are 0, end the simulation.
      if (H_t == 0 & V_t == 0) {
        H <- append(H, 0)
        V <- append(V, 0)
        time <- append(time, max_time)
        break
      }
    }
    
    sim_row = data.frame(
      w = w, 
      R0 = R0, 
      sigma = sigma,
      time = time,
      H = H,
      V = V
    )
    sims_df = rbind(sims_df, sim_row)
  }
  sims_df
}


# Run simulations ----
sims_out = init_sims_df
for (q in (1:dim(variables_df)[1])) {
  print(q)
  # Basic Reproduction Number
  R0 <- variables_df$R0[q]
  # Environmental oise level
  sigma <- variables_df$sigma[q]
  
  # calculate Thv based on R0 value
  Thv <- (R0^2) / ((b^2 * Tvh * Nv) / (Nh * gammah * muv)) 
  
  
  sim_summary =  cbind(data.frame(w = as.integer()), init_dataframe)
  
  # set up vectors to store results with value at time=0
  H <- c(H0)
  V <- c(V0)
  time <- c(0)
  
  H_t <- H0
  V_t <- V0
  
  for (w in (1:num_sims)) {
    H_t <- H0
    V_t <- V0
    H <- c(H0)
    V <- c(V0)
    time = c(0)
    # calculate random values all at once to avoid repeated calls
    dW_df = tibble(dW1 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat)),
                   dW2 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat)),
                   dW3 = rnorm(num_timesteps, mean = 0, sd = sqrt(deltat))
    )
    
    for (time_index in 1:length(time_vec)) {
      # following equations in paper draft-double check these ****
      dW1 = dW_df$dW1[time_index]
      dW2 = dW_df$dW2[time_index]
      dW3 = dW_df$dW3[time_index]
      
      H_s <- H_t
      V_s <- V_t
      
      dH_det = ((b * Thv / Nh) * V_s * (Nh - H_s) - gammah * H_s)
      dH_demo = sqrt(((b * Thv / Nh) * V_s * (Nh - H_s) + gammah * H_s))
      dH_enviro = alphab * (Thv / Nh) * V_s * (Nh - H_s)
      
      dH <- dH_det * deltat + dH_demo * dW1 + sigma * dH_enviro * dW3
      H_t = max(H_s + dH, 0)
      
      dV_det = ((b * Thv / Nh) * V_s * (Nh - H_s) - gammah * H_s)
      dV_demo = sqrt(((b * Thv / Nh) * V_s * (Nh - H_s) + gammah * H_s))
      dV_enviro = alphab * (Thv / Nh) * V_s * (Nh - H_s)
      
      dV <- dV_det * deltat + dV_demo * dW2 + sigma * dV_enviro * dW3
      V_t = max(V_s + dV, 0)
      
      # Add the values to vector
      H <- append(H, H_t)
      V <- append(V, V_t)
      time <- append(time, time_vec[time_index])
      
      # if both H_t and V_t are 0, end the simulation.
      if (H_t == 0 & V_t == 0) {
        H <- append(H, 0)
        V <- append(V, 0)
        time <- append(time, max_time)
        break
      }
    }
    
    sim_row = data.frame(
      w = w, 
      R0 = R0, 
      sigma = sigma,
      time = time,
      H = H,
      V = V
    )
    
    sims_out = rbind(sims_out, sim_row)
  }
}


# Calculate summary statistics
summary_df = sims_out %>% 
  # Calculate statistics for each simulation
  summarise(
    .by = c(R0, sigma, w),
    endemic = H[time == max_time] > 0, #
    small_outbreak = any(H > 10),
    big_outbreak = any(H > 100),
    max_cases = max(H),
    duration = ifelse(
      max_cases > 1, # If an outbreak actually occurred...
      max(time[H > 1]) - min(time[H > 1]), # Measure how long cases were above 1
      0 # Otherwise the duration is zero
    )
  ) %>% 
  # Summarise for each combination of R0 and sigma
  pivot_longer(cols = endemic:duration) %>% 
  mutate(
    .by = c(R0, sigma, name),
    mean = mean(value),
    variance = var(value),
    sum_value = as.integer(sum(value)) # needed for Clopper-Pearon confidence intervals
  ) %>% 
  select(-c(w, value)) %>% distinct() %>%
  rowwise() %>% 
  mutate(
    lower_ci = ifelse(
      name %in% c("max_cases", "duration"),
      max(0, mean - 1.96 * sqrt(variance)),
      clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "lower")$Lower.limit
    ),
    upper_ci = ifelse(
      name %in% c("max_cases", "duration"),
      min(Nh, mean + 1.96 * sqrt(variance)),
      clopper.pearson.ci(sum_value, as.integer(num_sims), alpha = 0.05, CI = "lower")$Upper.limit
    )
  )


# Plots ----
summary_df$R0 = factor(summary_df$R0, levels = rev(R0s))
sims_df$R0 = factor(sims_df$R0, levels = rev(R0s))

R0_colors = (c4a("tableau.classic_orange_blue", 20))[c(1:10, 16:20)]

# [] !!! set up nice plotting to match previous figures

# Generic plot function
generic_plot_function <- function(output_name) {
  summary_df %>% 
    filter(name == output_name) %>% 
    ggplot(aes(x = sigma, y = mean, group = R0, color = as.factor(R0))) +
    # Mean
    geom_line() +
    # Confidence intervals
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(R0)), alpha = 0.1) + 
    scale_color_manual("R0", values = R0_colors) +
    scale_fill_manual("R0", values = R0_colors) +
    scale_x_continuous(TeX("Environmental noise strength [$\\sigma]"),
                       limits = c(0,NA),
                       expand = c(0,0)) +
    scale_y_continuous(output_name,
                       limits = c(0,NA),
                       expand = c(0,0)) +
    theme_half_open()
}


# All simulations by R0 and sigma
All_sims_plot <- sims_out %>% 
  mutate(.by = c(w, R0, sigma),
         # Color according to whether the outbreak persists
         endemic = H[time == max_time] > 0) %>% 
  ggplot(aes(x = time/365, y = H, group = w, color = endemic)) +
  geom_path(alpha = 0.5) +
  facet_grid(rows = vars(R0), cols = vars(sigma), scales = 'free') +
  scale_color_manual(values = c("#DF536B", "#61D04F"), breaks = c("TRUE", "FALSE")) +
  scale_x_continuous("Time (years)", limits = c(0, NA)) +
  scale_y_continuous("Infected humans") +
  guides(color = "none") +
  theme_half_open() 

# Pr(endemic)
Pr_end_plot <- summary_df %>% 
  filter(name == "endemic") %>% 
  ggplot(aes(x = sigma, y = mean, group = R0, color = as.factor(R0))) +
  # Mean
  geom_line() +
  # Confidence intervals
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(R0)), alpha = 0.2) + 
  theme_half_open()

# Peak cases
Peak_cases_plot <-  summary_df %>% 
  filter(name == "max_cases") %>% 
  ggplot(aes(x = sigma, y = mean, group = R0, color = as.factor(R0))) +
  # Mean
  geom_line() +
  # Confidence intervals
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(R0)), alpha = 0.2) + 
  theme_half_open()

# Small outbreak
Small_outbreak_plot <- summary_df %>% 
  filter(name == "small_outbreak") %>% 
  ggplot(aes(x = sigma, y = mean, group = R0, color = as.factor(R0))) +
  # Mean
  geom_line() +
  # Confidence intervals
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(R0)), alpha = 0.2) + 
  theme_half_open()


# Big outbreak
Big_outbreak_plot <-  summary_df %>% 
  filter(name == "big_outbreak") %>% 
  ggplot(aes(x = sigma, y = mean, group = R0, color = as.factor(R0))) +
  # Mean
  geom_line() +
  # Confidence intervals
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(R0)), alpha = 0.2) + 
  theme_half_open()

# Duration
Duration_plot <- summary_df %>% 
  filter(name == "duration") %>% 
  ggplot(aes(x = sigma, y = mean, group = R0, color = as.factor(R0))) +
  # Mean
  geom_line() +
  # Confidence intervals
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(R0)), alpha = 0.2) + 
  theme_half_open()


