# Re-writing the SDE trajectory plotter to improve performance and reliability
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

# Drift terms
function f!(du,u,p,t)
    #(H,V) = u
    #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
    du[1] = (p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2] - p[6] * u[1]
    #du[1] = (b * τₕᵥ * (Nₕ - H) * V / Nₕ) - γₕ * H
    du[2] = (p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] - p[7] * u[2]
    #du[2] = (b * τᵥₕ * (Nᵥ - V) * H / Nₕ) - μᵥ * V]
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
    du[1,3] = p[8] * (p[1] * p[9] * p[2] * (p[4] - u[1]) / p[4]) * u[2]
    #du[1,3] = σ * αᵦ * τₕᵥ * (Nₕ - H) * V / Nₕ
    du[2,3] = p[8] * ((p[1] * p[9] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + (p[3] * p[10] * p[1] * (p[5] - u[2]) / p[4]) * u[1] + p[7] * p[11] * u[2])
    # Alt formulation: environmental noise causes a *proportional* change in parameters (1 + N) * baseline
    # du[1,3] = @views (p[8] * p[9] * p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2]
    # du[2,3] = @views p[8] * ((p[9] * p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + (p[10] * p[3] * p[1] * (p[5] - u[2]) / p[4]) * u[1] + p[11] * p[7] * u[2])
    #du[2,3] = σ * (αᵦ * τᵥₕ * H * (Nᵥ - V) / Nₕ + αₜ * b * H * (Nᵥ - V) / Nₕ + αₘ * V)
    nothing
end

### Define callbacks
# Callback 1: if infected host pop. goes negative, replace its value with zero
function condition1(u, t, integrator)
    u[1] .< 0.0
end
function affect1!(integrator)
    integrator.u[1] = 0.0
end
cb1 = DiscreteCallback(condition1, affect1!, save_positions=(false, false))
# Callback 2: if infected vector pop. goes negative, replace its value with zero
function condition2(u, t, integrator)
    u[2] .< 0.0
end
function affect2!(integrator)
    integrator.u[2] = 0.0
end
cb2 = DiscreteCallback(condition2, affect2!, save_positions=(false, false))

# Callback 4: if infected host pop. exceeds carrying capacity (10000), replace its value with carrying capacity
function condition4(u, t, integrator)
    u[1] - p[4] .> 0.0
end
function affect4!(integrator)
    integrator.u[1] = p[4]
end
cb4 = DiscreteCallback(condition4, affect4!, save_positions=(false, false))
# Callback 5:  if infected vector pop. exceeds carrying capacity (100000), replace its value with carrying capacity
function condition5(u, t, integrator)
    u[2] - p[5] .> 0.0
end
function affect5!(integrator)
    integrator.u[2] = p[5]
end
cb5 = DiscreteCallback(condition5, affect5!, save_positions=(false, false))
# Combine callbacks
#This callback terminates the simulation if the outbreak dies out defined as less than 1 host and vector infected
function terminate_condition(u,t,integrator)
    u[2] < 1.0 && u[1] < 1.0
end
function terminate_affect!(integrator)
    integrator.u[1] = 0.0
terminate!(integrator)
end
terminate_cb = DiscreteCallback(terminate_condition,terminate_affect!, save_positions=(false, false))
cbs = CallbackSet(cb1, cb2, cb4, cb5, terminate_cb)

## Set up ensemble problem
function ensemble_setup(drift, diffusion, parameters, callbacks)
    u0 = [0.0, 10.0]
    tspan = [0.0, 11*365.0]
    noise_rate_prototype = [0.0 0.0 0.0; 0.0 0.0 0.0]
    ensemble_problem = EnsembleProblem(SDEProblem(drift, diffusion, u0, tspan, parameters, noise_rate_prototype = noise_rate_prototype, callback = callbacks))
    return(ensemble_problem)
end

## Innermost solver
function inner_solver(ensemble_problem, trajectory_count::Int64)
    save_points = 0:10:(11*365)
    solution = @suppress solve(ensemble_problem, SRA1(), EnsembleSplitThreads(), trajectories = trajectory_count, verbose = false, saveat = save_points, save_everystep = false)
    return(DataFrame(solution))
end

## Conjoin data
function collate_solutions(parameter_grid, number_trajectories::Int64, drift, diffusion, callbacks, parameters::Vector{Float64} )
    bigdf = DataFrame()
    for s in ProgressBar(eachindex(parameter_grid))
        parameters[8] = parameter_grid[s][1] # sigma value
        R0 = parameter_grid[s][2]
        parameters[2] = tHV_from_R0(p, R0) # Thv value
        solution = inner_solver(ensemble_setup(drift, diffusion, parameters, callbacks), number_trajectories)
        states = solution.u
        times = solution.t
        for i in 1:size(states)[1]
            temp = reduce(vcat, transpose(states[i]))
            temp2 = reduce(hcat, [temp, times[i]])
            trajectory_indexes = repeat([i], size(temp)[1])
            sigma_val = repeat([p[8]], size(temp)[1])
            R0_val = repeat([R0], size(temp)[1])
            temp3 = reduce(hcat, [temp2, trajectory_indexes, sigma_val, R0_val])
            temp4 = DataFrame(temp3, [:H, :V, :t, :trajectory, :sigma, :R0])
        
            append!(bigdf,temp4)
          end

    end
    return(bigdf)
end

# Define parameters
function tHV_from_R0(p, R0) 
    (b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
    tHV = (R0.^2) / ((b.^2 * τᵥₕ * Nᵥ) / (Nₕ * γₕ * μᵥ))
    return(tHV)
  end
const sigmas = collect(0.0:0.0125:0.3) 
const R0s = [0.5; 0.75; 0.95:0.05:1.25; 2.0:1.0:6.0]# [0.75, 0.95, 1.05, 1.25, 2, 4, 6.5]
p = [0.3, 0.02450000, 0.9, 10000.0, 100000.0, 1/3, 1/21, 0.0, 1.0, 1.0, 1.0]
# const p = [0.3, 0, 0.5, 10000, 100000, 0.1, 0.1, 0, 1.0, 1.0, 1.0]
const THVs = tHV_from_R0(p, R0s)
# Create grid to reduce number of for loops
const parm_grid = collect(Iterators.product(sigmas, R0s))

# Precompile functions
collate_solutions(parm_grid[1:2], 2, f!, g!, cbs, p)

# Run and save data
final_data = collate_solutions(reverse(parm_grid), 100, f!, g!, cbs, p)



# cd("$(homedir())/Documents/Github/dahlin-noisy-mbds/results") do
cd("F:/GitHub/dahlin-noisy-mbds/results") do
    CSV.write("illustrate_trajectories_proportionalSDE.csv", final_data, transform = (col,val) -> something(val, missing))
end