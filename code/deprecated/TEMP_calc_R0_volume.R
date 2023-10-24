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
muV_star = 0.5

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

# Function: Calculate volume of region where R0 > 1
Vol_calc <- function(sigma_in) {
  
  b_hi = b_star + 2 * sigma_in #* b_star
  b_low = b_star - 2 * sigma_in #* b_star
  mu_hi = muV_star + 2 * sigma_in #* muV_star
  mu_low = muV_star - 2 * sigma_in #* muV_star
  tau_hi = tauVH_star + 2 * sigma_in #* tauVH_star
  tau_low = tauVH_star - 2 * sigma_in #* tauVH_star
  
  term_1 <- b_hi * tau_hi * 4 * sigma_in
  
  term_2 <- - (4/3) * C^(-1) * tau_hi^(1/2) * (mu_hi^(3/2) - mu_low^(3/2))
  
  term_3 <- 4 * C^(-2) * b_hi^(-1) * sigma_in * muV_star
  
  out <- term_1 + term_2 + term_3
  
  # normalize by the volume of the whole parameter space
  normal_factor = (4 * sigma_in)^3 # = (b_hi - b_low) * (mu_hi - mu_low) * (tau_hi - tau_low)
  
  out <- out / normal_factor
  
  return(out)
}
  
# 3) Define axes of variation numerically ---------------------------------

# sigma - strength of environmental noise

# Determine bounds on sigma so that parameters are well defined

# b_star - 2 * sigma > 0
# muV_star - 2 * sigma > 0
# tauVH_star - 2 * sigma > 0
# tauVH_star + 2 * sigma < 1
# 
# sigma < b_star / 2
# sigma < muV_star / 2
# sigma < tauVH_star / 2
# sigma < (1 - tauVH_star) / 2

max_sigma = min(b_star / 2, muV_star / 2, tauVH_star / 2, (1 - tauVH_star) / 2)

sigma_vals = seq(0, max_sigma, length.out = 1000)

# 4) Calculate Volume(R0 > 1) across axes ---------------------------------

out <- tibble(sigma = sigma_vals) %>% 
  mutate(Vol = Vol_calc(sigma))


# 5) Save data and visualizations -----------------------------------------

out %>% ggplot(aes(x = sigma, y = Vol)) +
  geom_line() +
  scale_y_log10() +
  xlab(unname(TeX("Environmental noise strength $(\\sigma)$"))) +
  ylab(unname(TeX("Volume of parameter space where $R_0 > 1$"))) +
  theme_cowplot()
