## Kyle Dahlin, University of Georgia, kydahlin@gmail.com
## Started October 2022
##
## Title: Stochastic Ross-Macdonald model ######################################
##
## Project: Noise-induced phase transitions: noise in the Ross-Macdonald model
##
## Purpose: Code and functions for numerically solving the Ross-Macdonald model
##          with noise added
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Define accessory functions
##           3) Model equations
##           4) SDE solver
##           5) Analysis and visualization functions
##
##
## Inputs:  
##             
##
## Outputs: 
##

# 1) Set-up,load packages, get data, etc.#######################################

## Load Libraries---------------------------------------------------------------
library(tidyverse)


# 2) Define accessory functions ################################################


# 3) Model equations ###########################################################

## If we assume that temperature T_t is a random variable, then we have to be 
## careful with how derive the system of SDEs. I'm still looking into how to 
## actually do this. In this case, we're then assuming that there is 
## environmental stochasticity (temperature) driving changes in demographic 
## parameters (which we determine using empirical TPCs).

# Parameters that may be random variables---------------------------------------
# lambda_V - mosquito recruitment rate (actually want to assume depends on larval parameters)
# sigma_V  - mosquito biting rate
# beta_V   - mosquito susceptibility (prob. of infection given contact)
# mu_V     - mosquito mortality rate
# eta_V    - pathogen development rate within mosquito
# f        - mosquito fecundity
# rho_L    - larval mosquito development rate
# mu_L     - larval mosquito mortaltiy rate

# Remaining parameters----------------------------------------------------------
# C_L      - larval mosquito carrying capacity
# lambda_H - host recruitment rate
# beta_H   - host susceptbility
# mu_H     - host mortality rate
# gamma_H  - host recovery rate


# 4) SDE solver ################################################################

# 5) Analysis and visualization functions ######################################