## Kyle Dahlin, University of Georgia, kydahlin@gmail.com
## Started September 2022
##
## Title: Mosquito trait functions #############################################
##
## Project: Noise-induced phase transitions: noise in the Ross-Macdonald model
##
## Purpose: Functions for determining mosquito traits as parameters for the 
##          transmission model. Traits may vary by mosquito species, pathogen
##          of interest, and temperature.
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Define accessory functions
##           3) Build parameter table-generating function for mosquito-pathogen pairs
##
##
## Inputs:  (currently not implemented but code can be modified to take 
##           different values for temperature)
##        temp_seq: a sequence of temperature values at which to evaluate
##          mosquito-pathogen traits
##             
##
## Outputs: MosquitoThermalResponse.csv
##          a table of vector+pathogen traits with columns of 
##          Species (= mosquito species), Pathogen (= name of pathogen), 
##          Temperature (in Celsius) and all mosquito-pathogen traits: 
##          maximum biting frequency (sigma_V), fecundity (f),
##          probability of egg survival to adulthood (delta_L),
##          immature mosquito development rate (rho_L), and adult mosquito 
##          mortality rate (mu_V) 

# 1) Set-up,load packages, get data, etc.#######################################

  ## Load Libraries---------------------------------------------------------------
library(tidyverse)

## Get data---------------------------------------------------------------------

# Read in mosquito-specific traits
MosqThermParms <- read_csv("./data/MosquitoThermalParameters.csv",
                           col_types = cols())


# Define mosquito species names and codes:
sp_names <- c('Aedes aegypti', 'Aedes albopictus','Culex quinquefasciatus', 
              'Anopheles spp.')
sp_code <- c('AE','AL','CQ', 'AN')
names(sp_names) <- sp_code

# Define pathogen names
path_names <- c('DENV', 'ZIKV', 'WNV', 'P.fal.')

# Define mosquito-pathogen pairs
## The generated parameter table should include a set of parameters for 
## each pair defined here
MosqPathPairs <- tibble(Mosquito_species=c('Aedes aegypti','Aedes albopictus',
                                           'Aedes aegypti',
                                           'Culex quinquefasciatus',
                                           'Anopheles spp.'),
                        Pathogen=c('DENV','DENV','ZIKV','WNV','P.fal.'))


# 2) Define accessory functions ################################################

## Define functions used to describe thermal response curves--------------------
# Briere function
Briere <- function(q,Tmin,Tmax) {
  function(t) {pmax(q*t*(t-Tmin)*(Tmax-t)^(1/2),0,na.rm = TRUE)}
}

# Quadratic
Quadratic <- function(q,Tmin,Tmax) {
  function(t) {pmax(-q*(t-Tmin)*(t-Tmax),0,na.rm=TRUE)}
}

# Linear
Linear <- function(q,z) {
  function(t) {pmax(-q*t+z,0,na.RM=FALSE)}
}


## Assign appropriate thermal response curve to traits--------------------------
### This is done through the following function:
get.thermal.response <- function(SpeciesName, PathogenName, ParameterName) {
  
  df <- filter(MosqThermParms, Species == SpeciesName, Pathogen == PathogenName,
               Parameter == ParameterName)
  
  # If parameter is not actually dependent on the pathogen,
  # replace 'PathogenName' with NA
  if (dim(df)[1]==0) df <- filter(MosqThermParms, Species == SpeciesName, 
                                  is.na(Pathogen), Parameter == ParameterName)
  
  # If the table contains more than one result, use the most recent
  if (dim(df)[1]>1) df <- filter(df,Year==max(Year))
  
  parms <- select(df,q,Tmin,Tmax)
  function_type <- select(df,FunctionType)
  
  if (dim(function_type)[1]==0) {return(function(t){NA})} else {
    
    if (function_type == 'Briere') {
      return(Briere(parms$q,parms$Tmin,parms$Tmax))
    }
    if (function_type == 'Quadratic') {
      return(Quadratic(parms$q,parms$Tmin,parms$Tmax))
    }
    if (function_type == 'Linear') {
      return(Linear(parms$q,parms$Tmin))
    }
  }
  
}



# 3) Build parameter table-generating function for mosquito-pathogen pairs #####

## Set up temperature vector----------------------------------------------------
# Lowest temperature
Temperature_minimum = 10
# Highest temperature
Temperature_maximum = 40
# Step size
Temperature_resolution = .0005

Temperature_vector = seq(Temperature_minimum,Temperature_maximum,
                         by = Temperature_resolution)

## NB: Care needs to be taken here as not all mosquito species and pathogens 
## have equivalent trait data available:
## *  EFD is computed differently for Cx. quinquefasciatus: EFD = ER * pO
## *  muV is computed differently for Anopheles spp.: muV = -log(p)
## *  deltaL is computed differently for Cx. quinquefasciatus: deltaL = eV*pLA
## *  betaV is computed differently for Culex quinquefasciatus+WNV and Aedes aegypti+ZIKV: betaV = bc

## Create vector trait table----------------------------------------------------
# Create a table with columns of:
# Mosquito_species, Pathogen, Temperature and Model parameters:
#    sigmaV, muV, fecundity,  deltaL, rhoL, etaV, betaV, deltaV
out <- tibble(Mosquito_species="", Pathogen = "", Temperature = numeric(),
              sigmaV = numeric(), fecundity = numeric(), deltaL = numeric(),
              rhoL = numeric(), muV = numeric(),
              etaV = numeric(), betaV = numeric())#, deltaV = numeric())

# Obtain and assign trait values------------------------------------------------
for (ii in 1:dim(MosqPathPairs)[1]) 
{
  Mosquito_species = MosqPathPairs$Mosquito_species[ii]
  Pathogen_name = MosqPathPairs$Pathogen[ii]
  # for (Temperature in Temperature_vector) 
  # {
  
  ## Accessory function for evaluating trait values
  get.parameter = function(ParameterName) {
    as.double(get.thermal.response(Mosquito_species, Pathogen_name, 
                                   ParameterName)(Temperature_vector))
  }
  ### Evaluate trait functions and assign parameter values
  
  ## Maximum biting frequency
  # sigmaV = 'BitingRate'
  sigmaV = get.parameter('BitingRate')
  
  ## Adult mortality rate
  # muV = 1/'lf', except for Anopheles spp: muV = -log(p)
  if (Mosquito_species == 'Anopheles spp.') {muV = -log(get.parameter('p'))}
  if (Mosquito_species == 'Aedes aegypti') {
    if (Pathogen_name == 'ZIKV') {
      muV = 1/get.thermal.response('Aedes aegypti', 'ZIKV', 'lf')(Temperature_vector)
    } else {
      muV = 1/get.thermal.response('Aedes aegypti', 'DENV', 'lf')(Temperature_vector)
      }
  }
  if (Mosquito_species != 'Anopheles spp.' & Mosquito_species != 'Anopheles spp.') {
    muV = 1/get.thermal.response(Mosquito_species, Pathogen_name, 'lf')(Temperature_vector) # to account for Ae. aegypti having add'l mortality when infected
  }
  
  ## Fecundity: number of female eggs laid per egg-laying event 
  # fecundity = 0.5*EFD/muV where EFD = 'EFD' except for Culex quinquefasciatus, EFD = 'ER' * 'pO'
  if (Mosquito_species == 'Culex quinquefasciatus') {
    EFD = get.parameter('ER') * get.parameter('pO')
  } else {
    EFD = get.parameter('EFD')
  }
  
  ## NB: If gonotrophic period length (1/sigmaV) exceeds lifespan (1/muV), 
  ##     set fecundity to zero
  fecundity = ifelse(1 / sigmaV > 1 / muV, 0, 0.5 * EFD / sigmaV) 
  
  ## Probability of survival from egg to adulthood
  # deltaL = 'pEA' except for Culex quinquefasciatus, deltaL = 'EV' * 'pLA'
  if (Mosquito_species == "Culex quinquefasciatus") {
    deltaL = get.parameter('EV') * get.parameter('pLA')
  } else {
    deltaL = get.parameter('pEA')
  }
  
  ## Mosquito development rate
  # rhoL = 'MDR'
  rhoL = get.parameter('MDR')
  
  ## Extrinsic  incubation rate
  # etaV = 'PDR'
  etaV = get.parameter('PDR')
  
  ## Vector competence
  # betaV = 'b' * 'c' except for Culex quinquefasciatus+WNV and Aedes aegypti+ZIKV, betaV = 'bc'
  if (Pathogen_name == "DENV") {
    betaV = get.parameter('b') * get.parameter('c')
  } else  {
    betaV = get.parameter('bc')
  }
   
  ## Add rows onto table
  out <- add_row(out, Mosquito_species = Mosquito_species, 
                 Pathogen = Pathogen_name, Temperature = Temperature_vector, 
                 sigmaV = sigmaV, fecundity = fecundity,
                 deltaL = deltaL, rhoL = rhoL, muV = muV, 
                 etaV = etaV, betaV = betaV)#, deltaV = deltaV)
}

# Save data frame---------------------------------------------------------------
write_csv(out,"results/MosquitoThermalResponse.csv")