# Noisy RM SDE solutions
# cd("$(homedir())\\GitHub\\NoisyMBDs")
#cd("/Users/karinebey/Documents/GitHub/NoisyMBDs") # !!! change to work across systems later
using CSV
using DataFrames
using DifferentialEquations
using NaNStatistics
using ProgressBars
using Statistics
using StatsBase
using Suppressor
using TimerOutputs
# using PlotlyJS

# Set constant to keep track of timer
const timer_output = TimerOutput()

# base function for mosquito borne disease
function f(du,u,p,t)
  #(H,V) = u
  #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
  du[1] = @views (p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2] - p[6] * u[1]
  #du[1] = (b * τₕᵥ * (Nₕ - H) * V / Nₕ) - γₕ * H
  du[2] = @views (p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] - p[7] * u[2]
  #du[2] = (b * τᵥₕ * (Nᵥ - V) * H / Nₕ) - μᵥ * V]

  #u[3] is a fake variable to track the maximum number of cases during a simulation
  du[3] = 0
  #u[4] is a fake variable to track the time when the maximum number of cases occurs during a simulation 
  du[4] = 0
  return(du)
end
# stochasticity terms
function g(du,u,p,t)
  #(V,H) = u
  #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ, σ, αᵦ, αₜ, αₘ) = p
  du[1,1] = @views sqrt(max(0,(p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2] + p[6] * u[1]))
  #du[1,1] = sqrt(max(0,(b * τₕᵥ * (Nₕ - H) * V / Nₕ) + γₕ * H))
  du[2,1] = 0
  du[1,2] = 0
  du[2,2] = @views sqrt(max(0,(p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + p[7] * u[2]))
  #du[2,2] = sqrt(max(0,(b * τᵥₕ * (Nᵥ - V) * H / Nₕ) + μᵥ * V))
  du[1,3] = @views (p[8] * p[9] * p[2] * (p[4] - u[1]) / p[4]) * u[2]
  #du[1,3] = σ * αᵦ * τₕᵥ * (Nₕ - H) * V / Nₕ
  du[2,3] = @views p[8] * ((p[9] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + (p[10] * p[1] * (p[5] - u[2]) / p[4]) * u[1] + p[11] * u[2])
  # Alt formulation: environmental noise causes a *proportional* change in parameters (1 + N) * baseline
  # du[1,3] = @views (p[8] * p[9] * p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2]
  # du[2,3] = @views p[8] * ((p[9] * p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + (p[10] * p[3] * p[1] * (p[5] - u[2]) / p[4]) * u[1] + p[11] * p[7] * u[2])
  #du[2,3] = σ * (αᵦ * τᵥₕ * H * (Nᵥ - V) / Nₕ + αₜ * b * H * (Nᵥ - V) / Nₕ + αₘ * V)
  du[3,1] = 0
  du[3,2] = 0
  du[3,3] = 0
  du[4,1] = 0
  du[4,2] = 0
  du[4,3] = 0
  return(du)
end

# Function: Calculate tHV values from specified R0 values
function tHV_from_R0(p, R0) 
  (b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
  tHV = (R0.^2) / ((b.^2 * τᵥₕ * Nᵥ) / (Nₕ * γₕ * μᵥ))
  return(tHV)
end

#Parameters: biting rate, THV, TVH, NH, NV, recovery, mortality, sigma, alphab, alphat, alpham
p = [0.3, 0, 0.5, 10000, 100000, 0.1, 0.1, 0, 0.0, 1.0, 1.0]


#environmental noise strength parameter 
const sigmas=collect(0:0.02:0.4) # sequence of numbers ranging from 0 to 0.4, incrementing by 0.02 #[0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4]
#Rate of transmission from hosts to vectors-We varied this parameter to get different R0 values
#Code for calculating Tvh for a given R0 value is in the main R code
# TVHs=[0.01250000, 0.02005556, 0.02450000, 0.03472222, 0.08888889,  0.35555556, 0.93888889]#R0=0.75, 0.95, 1.05, 1.25, 2, 4, 6.5
const R0s = [0.75, 0.95, 1.05, 1.25, 2, 4, 6.5]
const THVs = tHV_from_R0(p, R0s)
#[0.02005556, 0.02450000, 0.03472222, 0.08888889, 0.35555556, 0.93888889] # These are for R0=0.95, 1.05, 1.25, 2, 4, 6.5
#[ 0.02450000, 0.03472222, 0.05000000, 0.08888889, 0.35555556] #These are for R0=1.05, 1.25,1.5, 2, 4
#[0.01250000, 0.02005556, 0.02450000, 0.35555556] These are for R0=0.75, 0.95, 1.05, 4
#Now varying R0 by TVH since that is not a parameter we want to play with for sensitivity analysis

#initial conditions->10 infected mosquitos to star
const u0 = [0.0, 10.0, 0, 0]
#10 year simulation with time step of 1 day
const tspan = [0.0, 10*365]
const noise_rate_prototype = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
const dt = 1

#callback 1 makes Host population value 0 if goes negative
function condition1(u, t, integrator)
  u[1]<0
end
function affect1!(integrator)
  integrator.u[1] = 0
end
cb1 = DiscreteCallback(condition1, affect1!)
#callback 2 makes vector population value 0 if goes negative
function condition2(u, t, integrator)
  u[2]<0
end
function affect2!(integrator)
  integrator.u[2] = 0
end
cb2 = DiscreteCallback(condition2, affect2!)
#Callback 3 racks the maximum number of hosts infected and saves this value
function condition3(u, t, integrator)
  u[3]<u[1]
end
function affect3!(integrator)
  integrator.u[3] = integrator.u[1]
  integrator.u[4] = integrator.t
end
cb3 = DiscreteCallback(condition3, affect3!)
##callback 4 makes Host population value 10000 if it exceeds this
function condition4(u, t, integrator)
    u[1]> p[4] # 10000
end
function affect4!(integrator)
    integrator.u[1] = p[4] #10000
    integrator.u[3] = p[4] #10000
end
cb4 = DiscreteCallback(condition4, affect4!)
#callback 5 makes Vector population value 100000 if it exceeds this
function condition5(u, t, integrator)
  u[2]> p[5] # 100000
end
function affect5!(integrator)
  integrator.u[2] = p[5] #100000
end
cb5 = DiscreteCallback(condition5, affect5!)
#This callback terminates the simulation if the outbreak dies out defined as less than 1 host and vector infected
function terminate_affect!(integrator)
  integrator.u[1] = 0
  terminate!(integrator)
end
function terminate_condition(u,t,integrator)
  u[2]<1 && u[1]<1
end
terminate_cb = DiscreteCallback(terminate_condition,terminate_affect!)
cbs = CallbackSet(cb1, cb2, cb3, cb4, cb5, terminate_cb)

# Function: reduce to batches
function reduction(u,batch,I)
  tmp = sum(cat(batch..., dims = 5), dims = 5)/length(I)
  length(u) == 0 && return tmp, false
  cat(u, tmp, dims = 5), false
end

#Output function saves only the value of the variabes at the final time step
function output_func(sol, i)
  last(DataFrame(sol)), false
end

# Create grid to reduce number of for loops
const parm_grid = collect(Iterators.product(sigmas, THVs))

# Function: Calculate SDE solution trajectories
function SDE_solve_func(parms, num_runs, num_trajectories)
  #this makes an empty dataframe to fill with the final data
  df_oprob = DataFrame([Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
  Float64[],Float64[],Float64[], Float64[], Float64[], Float64[],
  Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]], 
  ["Thv", "sigma","prob_end_mean", "prob_end_25", "prob_end_75", "prob_out10_mean","prob_out10_25", "prob_out10_75", "prob_out100_mean","prob_out100_25", "prob_out100_75","max_cases_all_mean", 
  "max_cases_all_25", "max_cases_all_75", "max_cases_out_mean", "max_cases_out_25", "max_cases_out_75", "num_cases_mean","num_cases_25", "num_cases_75", 
  "end_time_mean", "end_time_25", "end_time_75", "end_time_out_mean", "end_time_out_25", "end_time_out_75",
  "time_max_mean", "time_max_25", "time_max_75"])
    
  # Initialize saved csv
  # CSV.write("julia_mean_test.csv", df_oprob)
  
  #In these loops we vary the parameters sigma and Tvh (R0) and complete 1000 runs of 1000 iterations for each parameter combination
  #Each run is sampled as a distribution to get mean and percentiles

  #loop for varying sigma and TVH values
  for s in ProgressBar(eachindex(parms))
    p[8]=parms[s][1] # sigma value
    p[2]=parms[s][2] # Thv value
    #Sets up empty dataframe for temporary data for each set of runs
    dftemp2 = DataFrame([Float64[],Float64[], Float64[], Float64[],Float64[], Float64[], Float64[],Float64[], Float64[]], ["eprob", "oprob10","oprob100", "max_cases_all", "max_cases_out", "num_cases", "end_time", "end_time_out", "time_max"])

      #Loop for each run
      # for r in 1:num_runs
        total_runs = num_runs * num_trajectories
        #defines the SDE problem to be solved
        @timeit timer_output "define_SDE" begin 
          prob = SDEProblem(f, g, u0, tspan, p, noise_rate_prototype=noise_rate_prototype, callback = cbs)
          # Set up an ensemble problem-run 1000 iterations and save only the final time step values
          ensembleprob = EnsembleProblem(prob, output_func=output_func)
        end
        @timeit timer_output "solve_SDE" begin 
          sol = DataFrame
          sol = @suppress DataFrame(solve(ensembleprob, EM(), dt = dt, EnsembleSplitThreads(), trajectories = total_runs, batch_size = num_runs, save_everystep = false, verbose = false))#, isoutofdomain = (u,p,t) -> any(x -> x < 0, u))
        end

        for run_index = 1:num_runs        
          row_indices = 10*(run_index-1) .+ (1:num_trajectories)
          temp_df = sol[row_indices, [1,2,4,5]]
          temp_df = rename(temp_df, :timestamp => :end_time, :value1 => :num_cases, :value3 => :max_cases_all, :value4 => :time_max)

          #for each iteration determine whether the disease became endemic (defined as lasting the full 10 year time period)
          temp_df.eprob = temp_df.end_time
          temp_df.eprob .= ifelse.(temp_df.eprob .< 3650, 0, temp_df.eprob)
          temp_df.eprob .= ifelse.(temp_df.eprob .== 3650, 1, temp_df.eprob)

          #for each iteration determine whether an outbreak occured (defined as > 10 hosts infected at one time point)
          temp_df.oprob10 = temp_df.max_cases_all
          temp_df.oprob10 .= ifelse.(temp_df.oprob10 .< 11, 0, temp_df.oprob10)
          temp_df.oprob10 .= ifelse.(temp_df.oprob10 .> 10, 1, temp_df.oprob10)

          #for each iteration determine whether a large outbreak occured (defined as > 100 hosts infected at one time point)
          temp_df.oprob100 = temp_df.max_cases_all
          temp_df.oprob100 .= ifelse.(temp_df.oprob100 .< 101, 0, temp_df.oprob100)
          temp_df.oprob100 .= ifelse.(temp_df.oprob100 .> 100, 1, temp_df.oprob100)

          #if the disease did not become endemic, for each iteration, calculate the max cases and end time
          temp_df.max_cases_out = temp_df.max_cases_all
          temp_df.max_cases_out .= ifelse.(temp_df.end_time .>3649, NaN, temp_df.max_cases_out)
          temp_df.end_time_out = temp_df.end_time
          temp_df.end_time_out .= ifelse.(temp_df.end_time .>3649, NaN, temp_df.end_time_out)

          temp_df = combine(temp_df, All() .=> [nanmean], renamecols=false)

          dftemp2 = vcat(dftemp2, temp_df)

          temp_df = 0
          
        end

        # dftemp2 = vcat(dftemp2, (DataFrame([[mean(temp_df.eprob)],[mean(temp_df.oprob10)], [mean(temp_df.oprob100)],[mean(temp_df.max_cases)], [mean(temp_df.H)],[mean(temp_df.t)], [mean(temp_df.time_max)]], ["eprob", "oprob10","oprob100", "max_cases_all",  "num_cases", "end_time", "time_max"])))
      # end

      @timeit timer_output "saving solutions" begin
          #calculates the summary results of all runs for this parameter combination-mean and 25th & 75th percentiles for all variables tracked
          # Load the dataset up calculated for previous parameters
          # df_oprob = CSV.read("julia_mean_test.csv", DataFrame)
          # Add in the new data
          df_oprob = vcat(df_oprob, (DataFrame([[p[2]],[p[8]], [mean(dftemp2.eprob)], [percentile(dftemp2.eprob,25)], [percentile(dftemp2.eprob,75)], 
          [mean(dftemp2.oprob10)], [percentile(dftemp2.oprob10,25)], [percentile(dftemp2.oprob10,75)], [mean(dftemp2.oprob100)], [percentile(dftemp2.oprob100,25)],
          [percentile(dftemp2.oprob100,75)],[mean(dftemp2.max_cases_all)], [percentile(dftemp2.max_cases_all,25)], [percentile(dftemp2.max_cases_all,75)],  
          [nanmean(dftemp2.max_cases_out)], [nanpctile(dftemp2.max_cases_out, 25)], 
          [nanpctile(dftemp2.max_cases_out, 75)],
          [mean(dftemp2.num_cases)], [percentile(dftemp2.num_cases, 25)], [percentile(dftemp2.num_cases, 75)], 
          [mean(dftemp2.end_time)], [percentile(dftemp2.end_time, 25)], [percentile(dftemp2.end_time, 75)], [nanmean(dftemp2.end_time_out)], 
          [nanpctile(dftemp2.end_time_out, 25)], 
          [nanpctile(dftemp2.end_time_out, 75)], [mean(dftemp2.time_max)], 
          [percentile(dftemp2.time_max,25)], [percentile(dftemp2.time_max,75)]], 
              ["Thv", "sigma","prob_end_mean", "prob_end_25", "prob_end_75", "prob_out10_mean","prob_out10_25", "prob_out10_75", "prob_out100_mean","prob_out100_25", "prob_out100_75","max_cases_all_mean", 
              "max_cases_all_25", "max_cases_all_75", "max_cases_out_mean", "max_cases_out_25", "max_cases_out_75", "num_cases_mean", "num_cases_25", "num_cases_75", 
              "end_time_mean", "end_time_25", "end_time_75", "end_time_out_mean", "end_time_out_25", "end_time_out_75",
              "time_max_mean", "time_max_25", "time_max_75"])))
          dftemp2 = 0
          # Save data set
          # CSV.write("julia_mean_test.csv", df_oprob, transform = (col,val) -> something(val, missing))
          # Clear variable to free up memory, except for the last iteration
          # if s != lastindex(parm_grid)
          #   df_oprob = DataFrame([Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
          #   Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]], 
          #       ["Tvh", "sigma","prob_end_mean", "prob_end_25", "prob_end_75", "prob_out10_mean","prob_out10_25", "prob_out10_75", "prob_out100_mean","prob_out100_25", "prob_out100_75","max_cases_all_mean", 
          #       "max_cases_all_25", "max_cases_all_75", "num_cases_mean", "num_cases_25", "num_cases_75", "end_time_mean", "end_time_25", "end_time_75", 
          #       "time_max_mean", "time_max_25", "time_max_75"])
          # end
        end
    end
  #Returns the NaNs to NaNs
  return(df_oprob)
end


# Run simulations
parms = parm_grid
num_runs = 10
num_trajectories = 10

# Initialize saved csv
# df_oprob = DataFrame([Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
# Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]], 
#     ["Tvh", "sigma","prob_end_mean", "prob_end_25", "prob_end_75", "prob_out10_mean","prob_out10_25", "prob_out10_75", "prob_out100_mean","prob_out100_25", "prob_out100_75","max_cases_all_mean", 
#     "max_cases_all_25", "max_cases_all_75", "num_cases_mean", "num_cases_25", "num_cases_75", "end_time_mean", "end_time_25", "end_time_75", 
#     "time_max_mean", "time_max_25", "time_max_75"])

#init = false
#if init
#save_df = DataFrame([Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
#Float64[],Float64[],Float64[], Float64[], Float64[], Float64[],
# Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]], 
 #["Thv", "sigma","prob_end_mean", "prob_end_25", "prob_end_75", "prob_out10_mean","prob_out10_25", "prob_out10_75", "prob_out100_mean","prob_out100_25", "prob_out100_75","max_cases_all_mean", 
 #"max_cases_all_25", "max_cases_all_75",  "max_cases_out_mean", "max_cases_out_25", "max_cases_out_75", "num_cases_mean", "num_cases_25", "num_cases_75", 
 #"end_time_mean", "end_time_25", "end_time_75", "end_time_out_mean", "end_time_out_25", "end_time_out_75",
 #"time_max_mean", "time_max_25", "time_max_75"])
#CSV.write("julia_mean_test.csv", save_df)
#end


#small_jobs = (num_runs < 101 & num_trajectories < 101)
#if small_jobs
#  for i in 1:7
#    indices = (i-1)*21 .+ (1:21)
#    parms = parm_grid[indices]
#    df_oprob = SDE_solve_func(parms, num_runs, num_trajectories)#####
#
#    #save_df = CSV.read("julia_mean_test.csv", DataFrame)
#    save_df = vcat(save_df, df_oprob)
#  end
#end
#  CSV.write("julia_mean_test.csv", save_df)
#  print(["Jobs done. i  = ", i])
# for i in 1:7
#   indices = (i-1)*21 .+ (1:21)
#   parms = parm_grid[indices]
  df_oprob = SDE_solve_func(parm_grid, num_runs, num_trajectories)

#save_df = CSV.read("julia_mean_test.csv", DataFrame)
  # save_df = vcat(save_df, df_oprob)

# end
# CSV.write("julia_mean_test.csv", read_df)
#CSV.write("julia_mean_test.csv", save_df)
#print(["Jobs done. i  = ", i])


cd("/Users/karinebey/Documents/GitHub/dahlin-noisy-mbds/") do
 CSV.write("alpha_test_tm_on.csv", df_oprob, transform = (col,val) -> something(val, missing))
end

# cd("results") do
#   CSV.write("julia_mean_test.csv", df_oprob, transform = (col,val) -> something(val, missing))
# end

# using PlotlyJS
# plot(
#     df_oprob,
#     x=:sigma, y=:prob_end_mean, color=:Tvh, quartilemethod="exclusive", kind="line",
#     Layout(boxmode="group")
# )

# # plot(
# #     df_oprob,
# #     x=:sigma, y=:prob_out_mean, color=:Tvh, quartilemethod="exclusive", kind="line",
# #     Layout(boxmode="group")
# # )

# # plot(
# #     df_oprob,
# #     x=:sigma, y=:max_cases_mean, color=:Tvh, quartilemethod="exclusive", kind="line",
# #     Layout(boxmode="group")
# # )

# # plot(
# #     df_oprob,
# #     x=:sigma, y=:num_cases_mean, color=:Tvh, quartilemethod="exclusive", kind="line",
# #     Layout(boxmode="group")
# # )

# # plot(
# #     df_oprob,
# #     x=:sigma, y=:end_time_mean, color=:Tvh, quartilemethod="exclusive", kind="line",
# #     Layout(boxmode="group")
# # )
