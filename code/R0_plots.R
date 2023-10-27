## Kyle Dahlin, University of Georgia, kydahlin@gmail.com
## Started October 2023
##
## Title: Calculate volume of R0 > 1 ######################################
##
## Project: Noise-induced phase transitions: noise in the Ross-Macdonald model
##
## Purpose: Calculate the likelihood of R0 > 1 given a certain level of environmental noise
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Define parameters and functions
##           3) Define axes of variation numerically
##           4) Calculate Volume(R0 > 1) across axes
##           5) Save data and visualizations
##


# 1) Set-up,load packages, get data, etc. ---------------------------------

## Load Libraries
library(tidyverse)
library(cowplot)
library(latex2exp)

# 2) Define parameters and functions --------------------------------------
# Baseline values for parameters depending on environmental noise:
# b  - mosquito biting rate
# tau_VH   - mosquito susceptibility (prob. of infection given contact)
# mu_V     - mosquito mortality rate

b_star = 0.3
tauVH_star = 0.5
muV_star = 0.1

# Remaining parameters:
# N_V      - mosquito population density
# tau_HV   - host susceptibility
# gamma_H  - host recovery rate
# N_H      - host population density
NV = 100000
tauHV = 0.03472222 # Gives R0 = 1.25 at baseline  
gammaH = 0.1
NH = 10000

C = sqrt(tauHV * NV / (gammaH * NH))

baseline_R0 = C * b_star * sqrt(tauVH_star / muV_star)

# 3) Define axes of variation numerically ---------------------------------

# Samples from Weiner process
sample_res = 1000
sample_vec = tibble(sample = rnorm(sample_res), ID = seq(1:sample_res))

# sigma - strength of environmental noise
sigma_res = 1000
sigma_vec = seq(0, 0.4, length.out = sigma_res)

eps = .Machine$double.eps

# Put together data
data_table <- cross_join(tibble(sigma = sigma_vec), sample_vec) %>% 
  rowwise() %>% 
  mutate(C = sqrt(tauHV * NV / (gammaH * NH)),
         b = max(eps, b_star + sigma * sample), 
         tauVH = max(eps, tauVH_star + sigma * sample), 
         muV = max(1/60, muV_star + sigma * sample),
         R0 = b * C * sqrt(tauVH / muV)
  )

summary_table <- data_table %>% 
  group_by(sigma) %>% 
  summarise(mean_R0 = mean(R0),
            sd_R0 = sd(R0))

condition_table <- data_table %>% 
  group_by(sigma) %>% 
  summarise(prop_end = sum(R0 < 1) / n())

# 4) Calculate Volume(R0 > 1) across axes ---------------------------------



# 5) Save data and visualizations -----------------------------------------

# R0 as a random variate
summary_table %>%
  ggplot(aes(x = sigma, y = mean_R0)) +
  geom_line() +
  geom_ribbon(aes(ymin = pmax(0, mean_R0 - 1.96*sd_R0), 
                  ymax = mean_R0 + 1.96*sd_R0),
              alpha = 0.25) +
  geom_hline(yintercept = 1, color = "red") +
  scale_x_continuous(unname(TeX("Environmental noise strength $(\\sigma)$")),
                     expand = c(0,0)) +
  scale_y_continuous(unname(TeX("Basic reproduction number ($R_0$)")),
                     expand = c(0,0)) +
  theme_cowplot()
ggsave("./figures/R0vsnoise.png", height = 6, width = 8)


# Proportion of simulations with R0 > 1
condition_table %>%
  ggplot(aes(x = sigma, y = prop_end)) +
  geom_line() +
  scale_x_continuous(unname(TeX("Environmental noise strength $(\\sigma)$")),
                     expand = c(0,0)) +
  scale_y_continuous(unname(TeX("$Pr(R_0 > 1)$")),
                     expand = c(0,0),
                     limits = c(0,1)) +
  theme_cowplot()
ggsave("./figures/propR0greater.png", height = 6, width = 8)

