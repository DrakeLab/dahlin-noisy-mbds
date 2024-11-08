# Numerical simulation of a Ross-Macdonald reflected SDE

# Load libraries
using DifferentialEquations
using Plots
using DataFrames
using DifferentialEquations.EnsembleAnalysis
using Main.Threads
using CSV
using ProgressBars
# using DiffEqMonteCarlo

# [] Later set this up as a function to parameters can be varied

# Set simulation parameters 

const deltat = 1.0 # timestep for Euler-Maruyama (EM) algorithm
const maxtime = 3650.0 # final time point
const time_vec = 1:deltat:maxtime # vector of timepoints at which to evaluate SDE
const num_timesteps = length(time_vec) # total number of timesteps for EM algorithm
const timespan = (0.0, maxtime)
const noise_rate_prototype = [0.0 0.0 0.0; 0.0 0.0 0.0] # prototype for noise

# Define the parameters
const b = 0.3 # biting rate
const Tvh = 0.5 # to vector transmission probability
const Nh = 1E4 # total number of humans
const Nv = 1E5 # total number of vectors
const muv = 0.1 # vector mortality rate
const gammah = 0.1 # human recovery rate
const alphab = 1 # toggle: environmental stochasticity in b
const alphat = 1 # toggle: environmental stochasticity in Tvh
const alpham = 1 # toggle: environmental stochasticity in muv
const q = [b, Tvh, Nh, Nv, muv, gammah]

## Initial conditions
const H0 = 0 # initial infected humans
const V0 = 10 # initial infected mosquitoes
const u0 = [H0, V0]

# const R0s = [6, 5, 4, 3, 2, 1.25, 1.2, 1.15, 1.1, 1.05, 1, 0.95, 0.75, 0.5, 0] # values of R0 to consider
const R0s = 0:0.25:5

function Thv_from_R0(q, R0) 
    (b, Tvh, Nh, Nv, muv, gammah) = q
    Thv = (R0.^2) / ((b.^2 * Tvh * Nv) / (Nh * gammah * muv))
    return(Thv)
  end

Thvs = Thv_from_R0(q,R0s) # used to vary R0

const sigmas = 0:0.1:2 # levels of environmental noise

# Initialize dataframes, if necessary

# Define the Ross-Macdonald Reflected SDE 

## Deterministic part
function dF_det!(du,u,p,t)
    # Name the variables
    (H,V) = u
    (Thv, sigma) = p#, b, Tvh, Nh, Nv, muv, gammah, alphab, alphat, alpham) = p
    
    du[1] = @views (b * Thv * (Nh - H) / Nh) * V - gammah * H
    du[2] = @views b * (Tvh / Nh) * H * (Nv - V) - muv * V
    return(du)
  end

## Stochastic part
function dF_stoch!(du,u,p,t)
    # Name the variables
    (H,V) = u
    (Thv, sigma) = p #, b, Tvh, Nh, Nv, muv, gammah, alphab, alphat, alpham) = p
    # Demographic stochasticity
    du[1,1] = @views sqrt((b * Thv * (Nh - H) / Nh) * V + gammah * H)
    du[2,1] = 0
    du[1,2] = 0
    du[2,2] = @views sqrt(b * (Tvh / Nh) * H * (Nv - V) + muv * V)
    # Environmental stochasticity
    du[1,3] = @views sigma * alphab * b * (Thv / Nh) * V * (Nh - H)
    du[2,3] = @views sigma * alphab * b * (Tvh / Nh) * H * (Nv - V) + sigma * alphat * b  * (Tvh / Nh) * H * (Nv - V) + sigma * alpham * muv *V
    return(du)
  end


## Reflections
# --- In Julia these are done via callbacks
condition_H(u,t,integrator) = true
function affect_H!(integrator)
    if integrator.u[1] < 0.0f0
        integrator.u[1] = 0.0f0
    elseif integrator.u[1] > Nh
        integrator.u[1] = Nh
    end
end
cb_H = DiscreteCallback(condition_H, affect_H!; save_positions = (false, false))

### Reflect if H goes below zero
function condition_H0(u, t, integrator) # Event when condition(u,t,integrator) == 0
    u[1]
end

function affect_H0!(integrator)
    integrator.u[1] = 0
end
cb_H0 = ContinuousCallback(condition_H0, affect_H0!; save_positions = (false, true))

### Reflect if H goes above Nh
function condition_HNh(u, t, integrator)
    u[1] > Nh
end
function affect_HNh!(integrator)
    integrator.u[1] = Nh
end
cb_HNh = ContinuousCallback(condition_HNh, affect_HNh!; save_positions = (false, true))


condition_V(u,t,integrator) = true
function affect_V!(integrator)
    if integrator.u[2] < 0.0f0
        integrator.u[2] = 0.0f0
    elseif integrator.u[2] > Nv
        integrator.u[2] = Nv
    end
end
cb_V = DiscreteCallback(condition_V, affect_V!; save_positions = (false, false))

### Reflect if V goes below zero
function condition_V0(u, t, integrator) # Event when condition(u,t,integrator) == 0
    u[2]
end

function affect_V0!(integrator)
    integrator.u[2] = 0.0f0
end
cb_V0 = ContinuousCallback(condition_V0, affect_V0!; save_positions = (false, true))

### Reflect if V goes above Nv
function condition_VNv(u, t, integrator) # Event when condition(u,t,integrator) == 0
    u[2] > Nv
end
function affect_VNv!(integrator)
    integrator.u[2] = Nv
end
cb_VNv = ContinuousCallback(condition_VNv, affect_VNv!; save_positions = (false, true))

### Stop simulations if case counts are low enough
function condition_terminate(u, t, integrator)
    u[2] < 1.0f0 && u[1] < 1.0f0
end
affect_terminate!(integrator) = terminate!(integrator)
cb_terminate = ContinuousCallback(condition_terminate, affect_terminate!)


## Collect callbacks
# cbs = CallbackSet(cb_H0, cb_HNh, cb_V0, cb_VNv, cb_terminate)
cbs = CallbackSet(cb_H, cb_V, cb_terminate)

# Define parameter values to iterate over
parameter_values = [(Thv, sigma) for Thv in Thvs, sigma in sigmas]


# Get example trajectories for grid plot ----------------------------------------------------------------------------
# Initialize an empty DataFrame to store results
param_names = ["H", "V"]

function run_sims(num_runs, parameter_values)
    results = DataFrame(time = Float64[], Thv = Float64[], sigma = Float64[], run = Int[], H = Float64[], V = Float64[])


    for i in ProgressBar(eachindex(parameter_values))

        ## Set up ensemble SDE problem
        prob = SDEProblem(dF_det!, dF_stoch!, u0, timespan, parameter_values[i], noise_rate_prototype = noise_rate_prototype, callback = cbs)
        ensembleprob = EnsembleProblem(prob)
        ## Run SDE solver
        sol = solve(ensembleprob, EM(), dt = deltat, EnsembleThreads(); trajectories = num_runs)

        # Collect results into a tidy DataFrame
        for run_id in 1:num_runs
            trajectory = sol[run_id]  # Get the individual trajectory solution
            
            results_append = DataFrame(
                time = trajectory.t,
                Thv = fill(parameter_values[i][1], length(trajectory.t)),
                sigma = fill(parameter_values[i][2], length(trajectory.t)),
                run = fill(run_id, length(trajectory.t)),
                H = [u[1] for u in trajectory.u],
                V = [u[2] for u in trajectory.u]
            )
            
            # Append to the main results DataFrame
            append!(results, results_append)
        end

    end
    return results
end

# Collect summarized outputs ----------------------------------------------------------------------------
prob = SDEProblem(dF_det!, dF_stoch!, u0, timespan, parameter_values[1], noise_rate_prototype = noise_rate_prototype, callback = cbs)
ensembleprob = EnsembleProblem(prob)
## Run SDE solver
sol = solve(ensembleprob, EM(), dt = deltat, EnsembleThreads(); trajectories = 1)

function collect_outputs(num_runs, parameter_values)

    # Initialize a DataFrame to store results for each trajectory and parameter combination
    results = DataFrame(Thv = Float32[], sigma = Float32[], run = Int[], max_value = Float32[], exceeded_10 = Bool[],
    exceeded_100 = Bool[], positive_at_final = Bool[], positive_duration = Float32[])

    for i in ProgressBar(eachindex(parameter_values))

        # Get trajectories
        prob = SDEProblem(dF_det!, dF_stoch!, u0, timespan, parameter_values[i], noise_rate_prototype = noise_rate_prototype, callback = cbs)
        ensembleprob = EnsembleProblem(prob)
        ## Run SDE solver
        sol = solve(ensembleprob, EM(), dt = deltat, EnsembleThreads(); trajectories = num_runs)

        # Analyze each trajectory
        for run_id in 1:num_runs
            trajectory = sol[run_id]  # Get the individual trajectory solution

            # Extract the values and times for the second state variable (u[2])
            H_values = [u[1] for u in trajectory.u]
            V_values = [u[2] for u in trajectory.u]
            times = trajectory.t
            
            # Calculate the required statistics
            max_value = maximum(H_values)
            exceeded_10 = any(v > 10.0f0 for v in H_values)
            exceeded_100 = any(v > 100.0f0 for v in H_values)
            positive_at_final = H_values[end] > 1.0f0
            positive_duration = sum((t2 - t1) for (v1, v2, t1, t2) in zip(V_values[1:end-1], V_values[2:end], times[1:end-1], times[2:end]) if v1 > 0.0f0 || v2 > 0.0f0)
            Thv = parameter_values[i][1]
            sigma = parameter_values[i][2]

            # Append results for this trajectory to the DataFrame
            push!(results, (Thv, sigma, run_id, max_value, exceeded_10, exceeded_100, positive_at_final, positive_duration))
        end
    end
    return(results)
end

test = collect_outputs(1, parameter_values)

num_sims = 10000::Int# number of simulations to Run

collect_all = collect_outputs(num_sims, parameter_values)
CSV.write(joinpath(dirname(dirname(pwd())), "data", "collect_all_outputs.csv"), collect_all)

trajectories_for_grid_plot = run_sims(50, parameter_values)
CSV.write(joinpath(dirname(dirname(pwd())), "data", "trajectories_for_grid_plot.csv"), trajectories_for_grid_plot)