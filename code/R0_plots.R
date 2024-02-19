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
library(ggh4x)


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
tauHV_125 = 0.03472222
tauHV_65 = tauHV_125 * 6.5 / 1.25
tauHV_vals = c(tauHV_125, tauHV_65) # Gives R0 = 1.25 at baseline
gammaH = 0.1
NH = 10000

# C = sqrt(tauHV * NV / (gammaH * NH))

# baseline_R0 = C * b_star * sqrt(tauVH_star / muV_star)

# 3) Define axes of variation numerically ---------------------------------

# Samples from Weiner process
sample_res = 2 * 1000
sample_vec1 = tibble(sample = rnorm(sample_res/2), ID = seq(1:(sample_res/2)))
sample_vec2 = tibble(sample = rnorm(sample_res/2), ID = seq(1:(sample_res/2)))

# sigma - strength of environmental noise
sigma_res = 100
sigma_vec = c(seq(0, 0.4, length.out = sigma_res), 0, 0.1, 0.2, 0.3, 0.4) %>% 
  sort() %>% unique()

eps = .Machine$double.eps

# Put together data
data_table <- rbind(tibble(tauHV = tauHV_125, sample_vec1),
                    tibble(tauHV = tauHV_65, sample_vec2)) %>% 
  cross_join(tibble(sigma = sigma_vec)) %>%
  rowwise() %>%
  mutate(
    R0 = max(eps, b_star + sigma * sample) * sqrt(tauHV * NV / (gammaH * NH)) * sqrt( max(eps, tauVH_star + sigma * sample) / (max(1/60, muV_star + sigma * sample)))
    # C = sqrt(tauHV * NV / (gammaH * NH)),
         # b = max(eps, b_star + sigma * sample), 
         # tauVH = max(eps, tauVH_star + sigma * sample), 
         # muV = max(1/60, muV_star + sigma * sample),
         # R0 = b * C * sqrt(tauVH / muV)
  ) %>% 
  mutate(base_R0_label = case_when(
    tauHV == tauHV_125 ~ 1.25,
    tauHV == tauHV_65 ~ 6.5
  ))

summary_table <- data_table %>% 
  group_by(sigma, base_R0_label) %>% 
  summarise(mean_R0 = mean(R0),
            sd_R0 = sd(R0))

condition_table <- data_table %>% 
  group_by(sigma, base_R0_label) %>% 
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
  facet_wrap(~ base_R0_label, scales = "free", ncol = 1,
             labeller = as_labeller(function(string) {paste0("R0 = ", string)}))  +
  theme_cowplot() +
  theme(
    strip.background = element_rect(color = NA, fill = NA)
  )
ggsave("./figures/R0vsnoise.png", height = 6, width = 8)

# Discrete version
data_table %>%
  filter(sigma %in% c(0, 0.1, 0.2, 0.3, 0.4)) %>%
  mutate(sigma = as.factor(sigma)) %>%
  ggplot(aes(x = sigma, y = R0)) +
  geom_violin() +
  geom_line() +
  geom_hline(yintercept = 1, color = "red") +
  scale_x_discrete(unname(TeX("Environmental noise strength $(\\sigma)$")),
                   expand = c(0,0)) +
  scale_y_continuous(unname(TeX("Basic reproduction number ($R_0$)")),
                     expand = c(0,0)) +
  theme_cowplot()
ggsave("./figures/violin_R0vsnoise.png", height = 6, width = 8)

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

