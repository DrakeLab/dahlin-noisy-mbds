# Ross-Macdonald SDE sensitivity analysis
using DifferentialEquations
using Plots
#using PlotlyJS
using CSV
using DataFrames
using Statistics
using StatsBase
using LatinHypercubeSampling
using Suppressor # to suppress warnings from SDE solver
using ProgressBars # to gauge how long the solver will take
using TimerOutputs


# 1) Define any necessary functions

#For sensitivity analysis, I will be tracking:
#Prob endemic,prob outbreak, peak cases, and outbreak duration
#So I will not be tracking u[4] which was time of peak cases in the function to hopefully simplify things

#basic function
function f(du,u,p,t)

    #(V,H) = u
    #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
    du[1] = (p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2] - p[6] * u[1]
    #du[1] = (b * τₕᵥ * (Nₕ - H) * V / Nₕ) - γₕ * H
    du[2] = (p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] - p[7] * u[2]
    #du[2] = (b * τᵥₕ * (Nᵥ - V) * H / Nₕ) - μᵥ * V]
    du[3] = 0
end
  # Diffusion terms
function g(du,u,p,t)
    #(V,H) = u
    #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ, σ, αᵦ, αₜ, αₘ) = p
    du[1,1] = sqrt(max(0,(p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2] + p[6] * u[1]))
    #du[1,1] = sqrt(max(0,(b * τₕᵥ * (Nₕ - H) * V / Nₕ) + γₕ * H))
    du[2,1] = 0
    du[1,2] = 0
    du[2,2] = sqrt(max(0,(p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + p[7] * u[2]))
    #du[2,2] = sqrt(max(0,(b * τᵥₕ * (Nᵥ - V) * H / Nₕ) + μᵥ * V))
    du[1,3] = (p[8] * p[9] * p[2] * (p[4] - u[1]) / p[4]) * u[2]
    #du[1,3] = σ * αᵦ * τₕᵥ * (Nₕ - H) * V / Nₕ
    du[2,3] = p[8] * ((p[9] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + (p[10] * p[1] * (p[5] - u[2]) / p[4]) * u[1] + p[11] * u[2])
    #du[2,3] = σ * (αᵦ * τᵥₕ * H * (Nᵥ - V) / Nₕ + αₜ * b * H * (Nᵥ - V) / Nₕ + αₘ * V)
    du[3,1] = 0
    du[3,2] = 0
    du[3,3] = 0
  
end

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


#I will randomly choose R0=1.25 for Thv value
#I added a fake variable to the end to use for determining the significance of prcc values
p = [0, 0.03472222, 0, 10000, 100000, 0.1, 0, 0.0, 1.0, 1.0, 1.0, 0]
const sigmas=collect(0:0.02:0.4)
const u0 = [0.0, 10.0, 0]
#10 year simulation with time step of 1 day
const tspan = [0.0, 10*365]
const noise_rate_prototype = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
const dt = 1
const timer_output = TimerOutput()

#driving function
function SDE_PRCC_func(par_list, num_runs, num_trajectories)
    #this makes an empty dataframe to fill with the final data
    
    df_prcc_results = DataFrame([Float64[], Float64[],Float64[],Float64[],Float64[], Float64[], Float64[], Float64[],Float64[], Float64[], Float64[],Float64[],Float64[],
    Float64[], Float64[],Float64[],Float64[], Float64[], Float64[],Float64[],Float64[]], 
  ["sigma", "PRCC_b_o10prob", "PRCC_muv_o10prob", "PRCC_tvh_o10prob", "PRCC_fake_o10prob", "PRCC_b_o100prob", "PRCC_muv_o100prob", "PRCC_tvh_o100prob",
  "PRCC_fake_o100prob","PRCC_b_eprob", "PRCC_muv_eprob", "PRCC_tvh_eprob", "PRCC_fake_eprob", "PRCC_b_peak", "PRCC_muv_peak", "PRCC_tvh_peak", "PRCC_fake_peak",
  "PRCC_b_dur", "PRCC_muv_dur", "PRCC_tvh_dur", "PRCC_fake_dur"])
      
  #Run over the 21 different sigma values
    for q in 1:21
#initialize a temporary data frame for each sigma value run
      df_oprob = DataFrame([Float64[], Float64[], Float64[],Float64[],Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]], 
    ["sigma", "b", "muv", "tvh", "fake", "eprob", "o10prob", "o100prob","peak_cases", "duration"])

      p[8] = sigmas[q]

    
    #In these loops we vary the parameters sigma and Tvh (R0) and complete 1000 runs of 1000 iterations for each parameter combination
    #Each run is sampled as a distribution to get mean and percentiles
  
    #loop for varying sigma and TVH values
      for s in ProgressBar(1:Int64(length(par_list)/4))
        p[1] = par_list[s] #biting rate
        p[3] = par_list[s+200] #Tvh
        p[7] = par_list[s+100] #mosquito mortality
        p[12] = par_list[s+300] #fake variable

        
        
          total_runs = num_runs * num_trajectories
          #defines the SDE problem to be solved
          @timeit timer_output "define_SDE" begin 
            prob = SDEProblem(f, g, u0, tspan, p, noise_rate_prototype=noise_rate_prototype, callback = cbs, reduction = reduction)
            # Set up an ensemble problem-run 1000 iterations and save only the final time step values
            ensembleprob = EnsembleProblem(prob, output_func=output_func) 
          end
          @timeit timer_output "solve_SDE" begin 
            sol = DataFrame
            sol = @suppress DataFrame(solve(ensembleprob, EM(), dt = dt, EnsembleSplitThreads(), trajectories = total_runs, batch_size = num_runs, save_everystep = false, verbose = false))#, isoutofdomain = (u,p,t) -> any(x -> x < 0, u))
          end
  
        temp_df = sol[!, [1,4]]
        temp_df = rename(temp_df, :timestamp => :duration, :value3 => :peak_cases)
  
            #for each iteration determine whether the disease became endemic (defined as lasting the full 10 year time period)
            temp_df.eprob = temp_df.duration
            temp_df.eprob .= ifelse.(temp_df.eprob .< 3650, 0, temp_df.eprob)
            temp_df.eprob .= ifelse.(temp_df.eprob .== 3650, 1, temp_df.eprob)
  
            #for each iteration determine whether an outbreak occured (defined as > 10 hosts infected at one time point)
            temp_df.o10prob = temp_df.peak_cases
            temp_df.o10prob .= ifelse.(temp_df.o10prob .< 11, 0, temp_df.o10prob)
            temp_df.o10prob .= ifelse.(temp_df.o10prob .> 10, 1, temp_df.o10prob)

            temp_df.o100prob = temp_df.peak_cases
            temp_df.o100prob .= ifelse.(temp_df.o100prob .< 11, 0, temp_df.o100prob)
            temp_df.o100prob .= ifelse.(temp_df.o100prob .> 10, 1, temp_df.o100prob)

            temp_df[!, :b] .= p[1]
            temp_df[!, :muv] .= p[7]
            temp_df[!, :tvh] .= p[3]
            temp_df[!, :sigma] .= p[8]
            temp_df[!, :fake] .=p[12]

            #take mean over all simulations
  
            temp_df = combine(temp_df, All() .=> [mean], renamecols=false)
  
            df_oprob = vcat(df_oprob, temp_df)
  
            temp_df = 0
        
      end

      #rank each value

      ranked_temp_df = combine(df_oprob, All() .=> ordinalrank; renamecols = false)
#calculate the prccs for each parameter for each of the two metrics
      df_prcc_results = vcat(df_prcc_results, (DataFrame([[p[8]],
    [partialcor(ranked_temp_df.b, ranked_temp_df.o10prob, Matrix(select(ranked_temp_df, [:muv, :tvh, :fake])))],
    [partialcor(ranked_temp_df.muv, ranked_temp_df.o10prob, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
    [partialcor(ranked_temp_df.tvh, ranked_temp_df.o10prob, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))], 
    [partialcor(ranked_temp_df.fake, ranked_temp_df.o10prob, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))],
    [partialcor(ranked_temp_df.b, ranked_temp_df.o100prob, Matrix(select(ranked_temp_df, [:muv, :tvh, :fake])))],
    [partialcor(ranked_temp_df.muv, ranked_temp_df.o100prob, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
    [partialcor(ranked_temp_df.tvh, ranked_temp_df.o100prob, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))], 
    [partialcor(ranked_temp_df.fake, ranked_temp_df.o100prob, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))],
    [partialcor(ranked_temp_df.b, ranked_temp_df.eprob, Matrix(select(ranked_temp_df, [:muv, :tvh,:fake])))],
    [partialcor(ranked_temp_df.muv, ranked_temp_df.eprob, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
    [partialcor(ranked_temp_df.tvh, ranked_temp_df.eprob, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))],
    [partialcor(ranked_temp_df.fake, ranked_temp_df.eprob, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))],
    [partialcor(ranked_temp_df.b, ranked_temp_df.peak_cases, Matrix(select(ranked_temp_df, [:muv, :tvh, :fake])))],
    [partialcor(ranked_temp_df.muv, ranked_temp_df.peak_cases, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
    [partialcor(ranked_temp_df.tvh, ranked_temp_df.peak_cases, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))], 
    [partialcor(ranked_temp_df.fake, ranked_temp_df.peak_cases, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))],
    [partialcor(ranked_temp_df.b, ranked_temp_df.duration, Matrix(select(ranked_temp_df, [:muv, :tvh, :fake])))],
    [partialcor(ranked_temp_df.muv, ranked_temp_df.duration, Matrix(select(ranked_temp_df, [:b, :tvh, :fake])))], 
    [partialcor(ranked_temp_df.tvh, ranked_temp_df.duration, Matrix(select(ranked_temp_df, [:b, :muv, :fake])))], 
    [partialcor(ranked_temp_df.fake, ranked_temp_df.duration, Matrix(select(ranked_temp_df, [:b, :muv, :tvh])))]], 
    ["sigma", "PRCC_b_o10prob", "PRCC_muv_o10prob", "PRCC_tvh_o10prob", "PRCC_fake_o10prob", "PRCC_b_o100prob", "PRCC_muv_o100prob", "PRCC_tvh_o100prob",
    "PRCC_fake_o100prob","PRCC_b_eprob", "PRCC_muv_eprob", "PRCC_tvh_eprob", "PRCC_fake_eprob", "PRCC_b_peak", "PRCC_muv_peak", "PRCC_tvh_peak", "PRCC_fake_peak",
    "PRCC_b_dur", "PRCC_muv_dur", "PRCC_tvh_dur", "PRCC_fake_dur"])))

      df_oprob = 0
      ranked_temp_df = 0

    end
    return(df_prcc_results)
end



# 2) Set up design matrix

# Set the distributions for each variable parameter
# KD: for now, let's use uniform distributions. We just need to decide on ranges. These are temporary, let's talk in more detail about these on Monday
const range_b = (1/14, 0.45) # biting rate ~= 1 / oviposition cycle length
const range_muv = (1/28, 2/10) # mortality rate ~= 1 / average lifespan
const range_tvh = (0.05, 1)
const range_fake = (1, 100) # vector competency is a probability. stay away from zero

# Latin Hypercube Sampling from the ranges described above

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

# 3) Calculate outputs

# We want a data frame with the following columns:
# b, muV, t_HV, outbreak probability, final number of infecteds, time to extinction
num_runs=1000
num_trajectories=100

df_prcc_results = SDE_PRCC_func(LHSamples, num_runs, num_trajectories)

plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_o10prob,df_prcc_results.PRCC_muv_o10prob, df_prcc_results.PRCC_tvh_o10prob, df_prcc_results.PRCC_fake_o10prob], label=["b" "muv" "tvh" "fake"])
plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_o100prob,df_prcc_results.PRCC_muv_o100prob, df_prcc_results.PRCC_tvh_o100prob, df_prcc_results.PRCC_fake_o100prob], label=["b" "muv" "tvh" "fake"])
plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_eprob,df_prcc_results.PRCC_muv_eprob, df_prcc_results.PRCC_tvh_eprob, df_prcc_results.PRCC_fake_eprob], label=["b" "muv" "tvh" "fake"])
plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_peak,df_prcc_results.PRCC_muv_peak, df_prcc_results.PRCC_tvh_peak, df_prcc_results.PRCC_fake_peak], label=["b" "muv" "tvh" "fake"])
plot(df_prcc_results.sigma, [df_prcc_results.PRCC_b_dur,df_prcc_results.PRCC_muv_dur, df_prcc_results.PRCC_tvh_dur, df_prcc_results.PRCC_fake_dur], label=["b" "muv" "tvh" "fake"])


cd("/Users/karinebey/Documents/GitHub/NoisyMBDs") do
  CSV.write("sensitivity_analysis_new_all.csv", df_prcc_results, transform = (col,val) -> something(val, missing))
end



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













