# Ross-Macdonald SDE sensitivity analysis
using Plots
using StatsBase
using LatinHypercubeSampling
using CodecZlib

include("SDE_sensitivity.jl")

function save_gzip(data, filename)
    open(joinpath(dirname(dirname(pwd())), "data", filename), "w") do io
        gzip_io = GzipCompressorStream(io)
        CSV.write(gzip_io, data, append = true, writeheader = true)
        close(gzip_io)
    end
end
# A sensitivity analysis of key model outputs to key transmission parameters

# Key parameters:
# - host population size (NH)
# - vector population size (NV)
# - transmission probability, vector to host (tau_VH)
# - biting rate (b)

# Key outputs:
# - probability of a large outbreak, Pr(outbreak > 100 cases in 10 years)
# - intensity of outbreak, peak number of cases attained in 10 years
# - duration of outbreak, final time when there are infections in at least one host and one vector

# Additional parameters of interest:
# - environmental noise intensity (sigma)

# Notes:
# * Keep the same number of simulations: 100,000
# * Must be done much more coarsely than the initial analysis
#   * 5 (or fewer) values of environmental noise: 0, 0.5, 1.0, 1.5, 2.0
#   * 5 (or fewer) values of R0 (via tau_HV): 0.95, 1.05, 2, 3, 4
#   * Take 5 (or fewer) values for each key parameter - 2 below baseline, 2 above
# * Consider plotting the distribution if feasible

# Coarse distributions of key outputs for varied parameters
#   For each noise level and R0 value, plot the distribution of key outputs for parameters varied one-at-a-time


## FUNCTIONS ##

## PARAMETER SETUP ##
## Set up parameter values and ranges

# Function: Get n values in increments to up to x% below and above the baseline value for each parameter
# Range is set to be 1-x% to 1+x% of the baseline value, with n values in between
function get_parameter_values(baseline, x, n)
    lower_bound = baseline * (1 - x)
    upper_bound = baseline * (1 + x)
    return range(lower_bound, upper_bound, length=n)
end


# Define values for each parameter to be used in the simulations
vary_percentage = 0.5 # vary parameters by 50% above and below baseline
const b_vals = get_parameter_values(b_default, vary_percentage, 5)
const Tvh_vals = get_parameter_values(Tvh_default, vary_percentage, 5)
const Nv_vals = get_parameter_values(Nv_default, vary_percentage, 5)
const Nh_vals = get_parameter_values(Nh_default, vary_percentage, 5)
const R0_vals = [0.95, 1.0, 1.05, 2.0, 3.0, 4.0, 5.0]
const sigma_vals = 0f0:0.1f0:2f0#[0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]

coarse_parameter_values = [(b_val, Tvh_val, Nh_val, Nv_val, Thv_from_R0(b_default, Tvh_default, Nh_default, Nv_default, R0), sigma) 
                    for b_val in b_vals, Tvh_val in Tvh_vals, Nh_val in Nh_vals, Nv_val in Nv_vals, R0 in R0_vals, sigma in sigma_vals]

## ANALYSIS: Output distributions ##
## Populate table with outputs for each parameter combination

## Calculate and save
print("Sensitivity analysis simulations")
collect_all_sensitivity = collect_outputs(dF_det_no_demo!, dF_stoch_no_demo!, 1_000, coarse_parameter_values; output_file=joinpath(dirname(dirname(pwd())), "data", "collect_outputs_sensitivity.csv"))

