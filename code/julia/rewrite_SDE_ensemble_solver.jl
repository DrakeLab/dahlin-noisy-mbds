# Re-writing the SDE ensemble solver code to improve performance and reliability
using CSV, DataFrames, DifferentialEquations, NaNStatistics, ProgressBars, Statistics, StatsBase, Suppressor, TimerOutputs

# Performance tips:
# 1. Performance critical code should be inside a function
#    - The functions should take arguments, instead of operating directly on global variables, see the next point.
#    - Variables should be local, or passed as arguments to functions, whenever possible.
# 2. Avoid untyped global variables
#    - global names are frequently constants, and declaring them as such greatly improves performance
# 3. Annotate values taken from untyped locations
#    - It is often convenient to work with data structures that may contain values of any type (arrays of type Array{Any}). But, if you're using one of these structures and happen to know the type of an element, it helps to share this knowledge with the compiler
#    - i.e. if you know the element of the array A should have type Int32, you write A[1]::Int32 in the function
#    - In the case that the type of A[1] is not known precisely, x can be declared via x = convert(Int32, A[1])::Int32 (Notice that convert itself needs a type declaration)
# 4. Break functions into multiple definitions
#    - Instead of having a function defined with different branches depending on the type of the argument, write function(argument::Type) = do_thing_to(argument) for each type of argument
# 5. Avoid changing the type of a variable
#    - If x is going to be type Float64 later, just initialize it with x::Float64 = 1, for example
# 6. Separate kernel functions (aka, function barriers)
#    - Many functions follow a pattern of performing some set-up work, and then running many iterations to perform a core computation. Where possible, it is a good idea to put these core computations in separate functions.
# 7. Access arrays in memory order, along columns
# 8. Pre-allocate outputs: If your function returns an Array or some other complex type, it may have to allocate memory. 
# 9. Fuse vectorized operations: Julia has a special dot syntax that converts any scalar function into a "vectorized" function call, and any operator into a "vectorized" operator, with the special property that nested "dot calls" are fusing: they are combined at the syntax level into a single loop, without allocating temporary arrays
# 10. Consider using views for slices: If you are doing many operations on the slice, this can be good for performance because it is more efficient to work with a smaller contiguous copy than it would be to index into the original array. On the other hand, if you are just doing a few simple operations on the slice, the cost of the allocation and copy operations can be substantial.
# 11. Consider StaticArrays.jl for small fixed-size vector/matrix operations
#     - If your application involves many small (< 100 element) arrays of fixed sizes (i.e. the size is known prior to execution), then you might want to consider using the StaticArrays.jl package. This package allows you to represent such arrays in a way that avoids unnecessary heap allocations and allows the compiler to specialize code for the size of the array, e.g. by completely unrolling vector operations (eliminating the loops) and storing elements in CPU registers.


# To create an array of size of size n with entries x, use fill(x, (n,n)). Make sure to declare the type of x.



# Define relevant functions

## Start with innermost solver and work outwards
function inner_solver(ensemble_problem, run_count::Int64)
        solution = solve(ensemble_problem, SRA1(), EnsembleSplitThreads(), trajectories = run_count, verbose = true, save_everystep = false)
        return(solution)
end

## Set up ensemble problem
function ensemble_setup(drift, diffusion, parameters, callbacks, output_function)
        u0 = [0.0, 10.0, 0, 0]
        tspan = [0.0, 10*365.0]
        noise_rate_prototype = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        prob = SDEProblem(drift, diffusion, u0, tspan, parameters, noise_rate_prototype = noise_rate_prototype, callback = callbacks)
        ensemble_problem = EnsembleProblem(prob, output_func = output_function)
        return(ensemble_problem)
end

## Collate solutions
function collate_solutions(solution)
        temp_df = solution[:, [1,2,4,5]]
        temp_df = rename(temp_df, :timestamp => :end_time, :value1 => :num_cases, :value3 => :max_cases_all, :value4 => :time_max)
        #for each iteration determine whether the disease became endemic (defined as lasting the full 10 year time period)
        temp_df.eprob .= ifelse.(temp_df.end_time .< 3650, 0, 1)
        #for each iteration determine whether an outbreak occured (defined as > 10 hosts infected at one time point)
        temp_df.oprob10 .= ifelse.(temp_df.max_cases_all .< 11, 0, 1)
        #for each iteration determine whether a large outbreak occured (defined as > 100 hosts infected at one time point)
        temp_df.oprob100 .= ifelse.(temp_df.max_cases_all .< 101, 0, 1)
        #if the disease did not become endemic, for each iteration, calculate the max cases and end time
        temp_df.max_cases_out .= ifelse.(temp_df.end_time .>3649, NaN, temp_df.max_cases_all)
        temp_df.end_time_out .= ifelse.(temp_df.end_time .>3649, NaN, temp_df.end_time)
        combine(temp_df, All() .=> [nanmean], renamecols=false)
end

# Main function
function main_function(parameter_grid, number_runs, number_trajectories, drift, diffusion, callbacks, output_function, parameters)
        # Preallocate output
        final_data = DataFrame(R0 = Float64[], sigma = Float64[],prob_end_mean = Float64[], prob_end_25 = Float64[], prob_end_75 = Float64[], prob_out10_mean = Float64[],prob_out10_25 = Float64[], prob_out10_75 = Float64[], prob_out100_mean = Float64[],prob_out100_25 = Float64[], prob_out100_75 = Float64[], max_cases_all_mean = Float64[], 
        max_cases_all_25 = Float64[], max_cases_all_75 = Float64[], max_cases_out_mean = Float64[], max_cases_out_25 = Float64[], max_cases_out_75 = Float64[], num_cases_mean = Float64[], num_cases_25 = Float64[], num_cases_75 = Float64[], 
        end_time_mean = Float64[], end_time_25 = Float64[], end_time_75 = Float64[], end_time_out_mean = Float64[], end_time_out_25 = Float64[], end_time_out_75 = Float64[],
        time_max_mean = Float64[], time_max_25 = Float64[], time_max_75 = Float64[])
        
        for s in ProgressBar(eachindex(parameter_grid))
                parameters[8] = parameter_grid[s][1] # sigma value
                parameters[2] = tHV_from_R0(p, parameter_grid[s][2]) # Thv value
                append!(final_data, batch_stats(parameter_grid[s][2], number_runs, number_trajectories, parameters, drift, diffusion, callbacks, output_function))
        end
        return(final_data)
end

function batch_stats(R0, number_runs, number_trajectories, parameters, drift, diffusion, callbacks, output_function)
        batch_data = DataFrame(Array{Union{Missing, Float64}}(missing, 0, 9), [:end_time, :num_cases, :max_cases_all, :time_max, :eprob, :oprob10, :oprob100, :max_cases_out, :end_time_out])

        for t in eachindex(1:number_runs)
                append!(batch_data, collate_solutions(DataFrame(inner_solver(ensemble_setup(drift, diffusion, parameters, callbacks, output_function),number_trajectories))))
        end
        
        DataFrame([[R0],[parameters[8]], [mean(batch_data.eprob)], [percentile(batch_data.eprob,25)], [percentile(batch_data.eprob,75)], 
          [mean(batch_data.oprob10)], [percentile(batch_data.oprob10,25)], [percentile(batch_data.oprob10,75)], [mean(batch_data.oprob100)], [percentile(batch_data.oprob100,25)],
          [percentile(batch_data.oprob100,75)],[mean(batch_data.max_cases_all)], [percentile(batch_data.max_cases_all,25)], [percentile(batch_data.max_cases_all,75)],  
          [nanmean(batch_data.max_cases_out)], [nanpctile(batch_data.max_cases_out, 25)], 
          [nanpctile(batch_data.max_cases_out, 75)],
          [mean(batch_data.num_cases)], [percentile(batch_data.num_cases, 25)], [percentile(batch_data.num_cases, 75)], 
          [mean(batch_data.end_time)], [percentile(batch_data.end_time, 25)], [percentile(batch_data.end_time, 75)], [nanmean(batch_data.end_time_out)], 
          [nanpctile(batch_data.end_time_out, 25)], 
          [nanpctile(batch_data.end_time_out, 75)], [mean(batch_data.time_max)], 
          [percentile(batch_data.time_max,25)], [percentile(batch_data.time_max,75)]], 
              ["R0", "sigma","prob_end_mean", "prob_end_25", "prob_end_75", "prob_out10_mean","prob_out10_25", "prob_out10_75", "prob_out100_mean","prob_out100_25", "prob_out100_75","max_cases_all_mean", 
              "max_cases_all_25", "max_cases_all_75", "max_cases_out_mean", "max_cases_out_25", "max_cases_out_75", "num_cases_mean", "num_cases_25", "num_cases_75", 
              "end_time_mean", "end_time_25", "end_time_75", "end_time_out_mean", "end_time_out_25", "end_time_out_75",
              "time_max_mean", "time_max_25", "time_max_75"])
end

### Output function saves only the value of the variables at the final time step
function output_func(sol, i)
        last(DataFrame(sol)), false
end

### Define callbacks
# Callback 1: if infected host pop. goes negative, replace its value with zero
function condition1(u, t, integrator)
        u[1] .< 0.0
end
function affect1!(integrator)
        integrator.u[1] = 0.0
end
cb1 = DiscreteCallback(condition1, affect1!)
# Callback 2: if infected vector pop. goes negative, replace its value with zero
function condition2(u, t, integrator)
        u[2] .< 0.0
end
function affect2!(integrator)
        integrator.u[2] = 0.0
end
cb2 = DiscreteCallback(condition2, affect2!)
#Callback 3 racks the maximum number of hosts infected and saves this value
function condition3(u, t, integrator)
        u[3] .< u[1]
end
function affect3!(integrator)
        integrator.u[3] = integrator.u[1]
        integrator.u[4] = integrator.t
end
cb3 = DiscreteCallback(condition3, affect3!)

# Callback 4: if infected host pop. exceeds carrying capacity (10000), replace its value with carrying capacity
function condition4(u, t, integrator)
        u[1] - p[4] .> 0.0
end
function affect4!(integrator)
        integrator.u[1] = p[4]
end
cb4 = DiscreteCallback(condition4, affect4!)
# Callback 5:  if infected vector pop. exceeds carrying capacity (100000), replace its value with carrying capacity
function condition5(u, t, integrator)
        u[2] - p[5] .> 0.0
end
function affect5!(integrator)
        integrator.u[2] = p[5]
end
cb5 = DiscreteCallback(condition5, affect5!)
# Combine callbacks
#This callback terminates the simulation if the outbreak dies out defined as less than 1 host and vector infected
function terminate_condition(u,t,integrator)
        u[2] < 1.0 && u[1] < 1.0
end
function terminate_affect!(integrator)
        integrator.u[1] = 0.0
terminate!(integrator)
end
terminate_cb = DiscreteCallback(terminate_condition,terminate_affect!)
cbs = CallbackSet(cb1, cb2, cb3, cb4, cb5, terminate_cb)

# Drift terms
function f!(du,u,p,t)
        #(H,V) = u
        #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
        du[1] = (p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2] - p[6] * u[1]
        #du[1] = (b * τₕᵥ * (Nₕ - H) * V / Nₕ) - γₕ * H
        du[2] = (p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] - p[7] * u[2]
        #du[2] = (b * τᵥₕ * (Nᵥ - V) * H / Nₕ) - μᵥ * V]

        #u[3] is a fake variable to track the maximum number of cases during a simulation
        du[3] = 0.0
        #u[4] is a fake variable to track the time when the maximum number of cases occurs during a simulation 
        du[4] = 0.0
        nothing
end
# Diffusion terms
function g!(du,u,p,t)
        #(H,V) = u
        #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ, σ, αᵦ, αₜ, αₘ) = p
        du[1,1] = sqrt((p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2] + p[6] * u[1])
        #du[1,1] = sqrt(max(0,(b * τₕᵥ * (Nₕ - H) * V / Nₕ) + γₕ * H))
        du[2,1] = 0.0
        du[1,2] = 0.0
        du[2,2] = sqrt((p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + p[7] * u[2])
        #du[2,2] = sqrt(max(0,(b * τᵥₕ * (Nᵥ - V) * H / Nₕ) + μᵥ * V))
        du[1,3] = (p[1] * p[8] * p[9] * p[2] * (p[4] - u[1]) / p[4]) * u[2]
        #du[1,3] = σ * αᵦ * τₕᵥ * (Nₕ - H) * V / Nₕ
        du[2,3] = p[8] * ((p[1] * p[9] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + (p[3] * p[10] * p[1] * (p[5] - u[2]) / p[4]) * u[1] + p[7] * p[11] * u[2])
        # Alt formulation: environmental noise causes a *proportional* change in parameters (1 + N) * baseline
        # du[1,3] = @views (p[8] * p[9] * p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2]
        # du[2,3] = @views p[8] * ((p[9] * p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + (p[10] * p[3] * p[1] * (p[5] - u[2]) / p[4]) * u[1] + p[11] * p[7] * u[2])
        #du[2,3] = σ * (αᵦ * τᵥₕ * H * (Nᵥ - V) / Nₕ + αₜ * b * H * (Nᵥ - V) / Nₕ + αₘ * V)
        du[3,1] = 0.0
        du[3,2] = 0.0
        du[3,3] = 0.0
        du[4,1] = 0.0
        du[4,2] = 0.0
        du[4,3] = 0.0
        nothing
end

# Define parameters
function tHV_from_R0(p, R0) 
        (b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
        tHV = (R0.^2) / ((b.^2 * τᵥₕ * Nᵥ) / (Nₕ * γₕ * μᵥ))
        return(tHV)
      end
const sigmas = collect(0.0:0.0125:0.3) 
const R0s = [0.5; 0.75; 0.95:0.05:1.25; 2.0:1.0:6.0]# [0.75, 0.95, 1.05, 1.25, 2, 4, 6.5]
# p = [0.3, 0.02450000, 0.9, 10000.0, 100000.0, 1/3, 1/21, 0.0, 1.0, 1.0, 1.0]
const p = [0.3, 0, 0.5, 10000, 100000, 0.1, 0.1, 0, 1.0, 1.0, 1.0]
const THVs = tHV_from_R0(p, R0s)
# Create grid to reduce number of for loops
const parm_grid = collect(Iterators.product(sigmas, R0s))
# Define function arguments

# Run functions to get data
# # Re-write to use R0 instead of Thv
final_data = main_function(parm_grid[1:2], 2, 2, f!, g!, cbs, output_func, p)
final_data = main_function(reverse(parm_grid), 1000, 100, f!, g!, cbs, output_func, p)

# Save data
cd("$(homedir())/Documents/Github/dahlin-noisy-mbds/results") do
        CSV.write("propSDE_means.csv", final_data, transform = (col,val) -> something(val, missing))
end