# Ross-Macdonald SDE sensitivity analysis
using Plots
using StatsBase
using LatinHypercubeSampling
using CodecZlib

include("mod_RM_Reflected_SDE.jl")

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

# Major analyses:
# 1) Coarse distributions of key outputs for varied parameters
#    a) For each noise level and R0 value, plot the distribution of key outputs for parameters varied one-at-a-time
# 2) PRCCs of key outputs to key parameters for each noise level and R0 value
#   a) Use a Latin Hypercube Sampling design to sample parameter space for each noise level and R0 value


## FUNCTIONS ##



## PARAMETER SETUP ##
## Set up parameter values and ranges

# Key parameters

# Noise values

# R0 values

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

## ANALYSIS 1: Output distributions ##
## Populate table with outputs for each parameter combination

## Calculate statistics

# 3) Calculate outputs

## Calculate and save
print("Sensitivity analysis simulations")
collect_all_sensitivity = collect_outputs(dF_det_no_demo!, dF_stoch_no_demo!, 1_000, coarse_parameter_values; output_file=joinpath(dirname(dirname(pwd())), "data", "collect_outputs_sensitivity.csv"))

## ANALYSIS 2: PRCCs ##

# Latin Hypercube Sampling from the ranges described above
## Set up design matrix of parameter values and outputs

# Set the ranges for each variable parameter
const range_b = (1/14, 0.45) # biting rate ~= 1 / oviposition cycle length, baseline = 0.3
# const range_muv = (1/28, 2/10) # mortality rate ~= 1 / average lifespan, baseline = 0.1
const range_Tvh = (0.05, 1) # baseline = 0.5
const range_Nv = (50_000, 150_000) # baseline = 100,000
const range_Nh = (5_000, 15_000) # baseline = 10,000
const range_fake = (1, 100) # dummy parameter range

# Run this to install necessary package
#using Pkg
#Pkg.add("LatinHypercubeSampling")


# Generate a random Latin Hypercube with d="number of parameter" dimensions and n="number of samples" sample points
const d = 4
const n = 100
LHS_plan = randomLHC(n, d)

# Scale the hypercube to ensure parameters fall in the ranges described above
LHSamples = scaleLHC(LHS_plan, [range_b, range_muv, range_tvh, range_fake])
LHSamples = vec(LHSamples)



# function SDE_PRCC_func(sigmas_indices, par_list, num_runs, num_trajectories)
#     #this makes an empty dataframe to fill with the final data
    
#     df_prcc_results = DataFrame([Float64[], Float64[],Float64[],Float64[],Float64[], Float64[], Float64[], Float64[],Float64[], Float64[], Float64[],Float64[],Float64[],
#     Float64[], Float64[],Float64[],Float64[], Float64[], Float64[],Float64[],Float64[]], 
#   ["sigma", "PRCC_b_o10prob", "PRCC_muv_o10prob", "PRCC_tvh_o10prob", "PRCC_fake_o10prob", "PRCC_b_o100prob", "PRCC_muv_o100prob", "PRCC_tvh_o100prob",
#   "PRCC_fake_o100prob","PRCC_b_eprob", "PRCC_muv_eprob", "PRCC_tvh_eprob", "PRCC_fake_eprob", "PRCC_b_peak", "PRCC_muv_peak", "PRCC_tvh_peak", "PRCC_fake_peak",
#   "PRCC_b_dur", "PRCC_muv_dur", "PRCC_tvh_dur", "PRCC_fake_dur"])
      
#   #Run over the 21 different sigma values
#     for q in sigmas_indices
#       print(q)
# #initialize a temporary data frame for each sigma value run
#       df_oprob = DataFrame([Float64[], Float64[], Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]], 
#     ["sigma", "b", "muv", "tvh", "fake", "eprob", "o10prob", "o100prob","peak_cases", "duration"])

#       p[8] = sigmas[q]

    
#     #In these loops we vary the parameters sigma and Tvh (R0) and complete 1000 runs of 1000 iterations for each parameter combination
#     #Each run is sampled as a distribution to get mean and percentiles
  
#     #loop for varying sigma and TVH values
#       for s in ProgressBar(1:Int64(length(par_list)/4))
#         p[1] = par_list[s] #biting rate
#         p[3] = par_list[s+200] #Tvh
#         p[7] = par_list[s+100] #mosquito mortality
#         p[12] = par_list[s+300] #fake variable

        
        
#           total_runs = num_runs * num_trajectories
#           #defines the SDE problem to be solved
#           @timeit timer_output "define_SDE" begin 
#             prob = SDEProblem(f, g, u0, tspan, p, noise_rate_prototype=noise_rate_prototype, callback = cbs, reduction = reduction)
#             # Set up an ensemble problem-run 1000 iterations and save only the final time step values
#             ensembleprob = EnsembleProblem(prob, output_func=output_func) 
#           end
#           @timeit timer_output "solve_SDE" begin 
#             sol = DataFrame
#             sol = @suppress DataFrame(solve(ensembleprob, EM(), dt = dt, EnsembleSplitThreads(), trajectories = total_runs, batch_size = num_runs, save_everystep = false, verbose = false))#, isoutofdomain = (u,p,t) -> any(x -> x < 0, u))
#           end
  
#         temp_df = sol[!, [1,4]]
#         temp_df = rename(temp_df, :timestamp => :duration, :value3 => :peak_cases)
  
#             #for each iteration determine whether the disease became endemic (defined as lasting the full 10 year time period)
#             temp_df.eprob = temp_df.duration
#             temp_df.eprob .= ifelse.(temp_df.eprob .< 3650, 0, temp_df.eprob)
#             temp_df.eprob .= ifelse.(temp_df.eprob .== 3650, 1, temp_df.eprob)
  
#             #for each iteration determine whether an outbreak occured (defined as > 10 hosts infected at one time point)
#             temp_df.o10prob = temp_df.peak_cases
#             temp_df.o10prob .= ifelse.(temp_df.o10prob .< 11, 0, temp_df.o10prob)
#             temp_df.o10prob .= ifelse.(temp_df.o10prob .> 10, 1, temp_df.o10prob)

#             temp_df.o100prob = temp_df.peak_cases
#             temp_df.o100prob .= ifelse.(temp_df.o100prob .< 11, 0, temp_df.o100prob)
#             temp_df.o100prob .= ifelse.(temp_df.o100prob .> 10, 1, temp_df.o100prob)

#             temp_df[!, :b] .= p[1]
#             temp_df[!, :muv] .= p[7]
#             temp_df[!, :tvh] .= p[3]
#             temp_df[!, :sigma] .= p[8]
#             temp_df[!, :fake] .=p[12]

#             #take mean over all simulations
  
#             temp_df = combine(temp_df, All() .=> [mean], renamecols=false)
  
#             df_oprob = vcat(df_oprob, temp_df)
  
#             temp_df = 0
        
#       end

#       #rank each value

#       ranked_temp_df = combine(df_oprob, All() .=> ordinalrank; renamecols = false)
# #calculate the prccs for each parameter for each of the two metrics
#       df_prcc_results = vcat(df_prcc_results, (DataFrame([[p[8]],
#     [partialcor(ranked_temp_df.b, ranked_temp_df.o10prob, Matrix(select(ranked_temp_df, [:muv, :tvh, :fake])))],
#     [partialcor(ranked_temp_df.muv, ranked_temp_df.o10prob, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
#     [partialcor(ranked_temp_df.tvh, ranked_temp_df.o10prob, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))], 
#     [partialcor(ranked_temp_df.fake, ranked_temp_df.o10prob, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))],
#     [partialcor(ranked_temp_df.b, ranked_temp_df.o100prob, Matrix(select(ranked_temp_df, [:muv, :tvh, :fake])))],
#     [partialcor(ranked_temp_df.muv, ranked_temp_df.o100prob, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
#     [partialcor(ranked_temp_df.tvh, ranked_temp_df.o100prob, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))], 
#     [partialcor(ranked_temp_df.fake, ranked_temp_df.o100prob, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))],
#     [partialcor(ranked_temp_df.b, ranked_temp_df.eprob, Matrix(select(ranked_temp_df, [:muv, :tvh,:fake])))],
#     [partialcor(ranked_temp_df.muv, ranked_temp_df.eprob, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
#     [partialcor(ranked_temp_df.tvh, ranked_temp_df.eprob, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))],
#     [partialcor(ranked_temp_df.fake, ranked_temp_df.eprob, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))],
#     [partialcor(ranked_temp_df.b, ranked_temp_df.peak_cases, Matrix(select(ranked_temp_df, [:muv, :tvh, :fake])))],
#     [partialcor(ranked_temp_df.muv, ranked_temp_df.peak_cases, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
#     [partialcor(ranked_temp_df.tvh, ranked_temp_df.peak_cases, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))], 
#     [partialcor(ranked_temp_df.fake, ranked_temp_df.peak_cases, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))],
#     [partialcor(ranked_temp_df.b, ranked_temp_df.duration, Matrix(select(ranked_temp_df, [:muv, :tvh, :fake])))],
#     [partialcor(ranked_temp_df.muv, ranked_temp_df.duration, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
#     [partialcor(ranked_temp_df.tvh, ranked_temp_df.duration, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))], 
#     [partialcor(ranked_temp_df.fake, ranked_temp_df.duration, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))]], 
#     ["sigma", "PRCC_b_o10prob", "PRCC_muv_o10prob", "PRCC_tvh_o10prob", "PRCC_fake_o10prob", "PRCC_b_o100prob", "PRCC_muv_o100prob", "PRCC_tvh_o100prob",
#     "PRCC_fake_o100prob","PRCC_b_eprob", "PRCC_muv_eprob", "PRCC_tvh_eprob", "PRCC_fake_eprob", "PRCC_b_peak", "PRCC_muv_peak", "PRCC_tvh_peak", "PRCC_fake_peak",
#     "PRCC_b_dur", "PRCC_muv_dur", "PRCC_tvh_dur", "PRCC_fake_dur"])))

#       df_oprob = 0
#       ranked_temp_df = 0

#     end
#     return(df_prcc_results)
# end

# We want a data frame with the following columns:
# b, muV, t_HV, outbreak probability, final number of infecteds, time to extinction
num_runs=1000
num_trajectories=100

df_prcc_results = SDE_PRCC_func(sigmas_indices, LHSamples, num_runs, num_trajectories)

# CSV.write("sensitivity_analysis_new_all.csv", df_prcc_results)

read_df = CSV.read("sensitivity_analysis_new_all.csv", DataFrame)
save_df = vcat(read_df, df_prcc_results)

# cd("/Users/karinebey/Documents/GitHub/NoisyMBDs") do
CSV.write("sensitivity_analysis_new_all.csv", save_df, transform = (col,val) -> something(val, missing))
# end

df_prcc_results = CSV.read("sensitivity_analysis_new_all.csv", DataFrame)

# plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_o10prob,df_prcc_results.PRCC_muv_o10prob, df_prcc_results.PRCC_tvh_o10prob, df_prcc_results.PRCC_fake_o10prob], label=["b" "muv" "tvh" "fake"])
# png("o10_SA.png")
# plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_o100prob,df_prcc_results.PRCC_muv_o100prob, df_prcc_results.PRCC_tvh_o100prob, df_prcc_results.PRCC_fake_o100prob], label=["b" "muv" "tvh" "fake"])
# png("o100_SA.png")
# plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_eprob,df_prcc_results.PRCC_muv_eprob, df_prcc_results.PRCC_tvh_eprob, df_prcc_results.PRCC_fake_eprob], label=["b" "muv" "tvh" "fake"])
# png("endemic_SA.png")
# plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_peak,df_prcc_results.PRCC_muv_peak, df_prcc_results.PRCC_tvh_peak, df_prcc_results.PRCC_fake_peak], label=["b" "muv" "tvh" "fake"])
# png("peak_SA.png")
# plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_dur,df_prcc_results.PRCC_muv_dur, df_prcc_results.PRCC_tvh_dur, df_prcc_results.PRCC_fake_dur], label=["b" "muv" "tvh" "fake"])
# png("duration_SA.png")




# A pretend matrix for now, but final should look like this
#output_df = DataFrame((b = LHSamples[:,1], muv = LHSamples[:,2], tvh = LHSamples[:,3], eprob=rand(100), oprob = rand(100), peak_cases =  rand(100), duration =  rand(100)))

# Might have to separate these out
#oprob_df = select(output_df, Not([:num_cases, :end_time]))
#num_cases_df = select(output_df, Not([:oprob, :end_time]))
#end_time_df = select(output_df, Not([:oprob, :num_cases]))


# 4) Calculate PRCCs of each output to each parameter

# We should end up with three tables (one for each output) with three entries each (one for each parameter)

# Rank transform input-output matrix

#using Pkg; Pkg.add("StatsBase")


#ranked_output_df = combine(df_oprob, All() .=> ordinalrank; renamecols = false)

# Calculate outbreak probability PRCCs
# KD: there's definitely a cleaner way to do this with split-apply-combine, but this works for now
#prcc_oprob_b = partialcor(ranked_output_df.b, ranked_output_df.oprob, Matrix(select(ranked_output_df, [:muv, :tvh])))
#prcc_oprob_muv = partialcor(ranked_output_df.muv, ranked_output_df.oprob, Matrix(select(ranked_output_df, [:b, :tvh])))
#prcc_oprob_tvh = partialcor(ranked_output_df.tvh, ranked_output_df.oprob, Matrix(select(ranked_output_df, [:b, :muv])))


#oprob_prccs = combine(ranked_output_df, [:b, :oprob] => ((parm,output) -> (prcc = partialcor(parm, output))))













