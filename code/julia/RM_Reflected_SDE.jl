# Numerical simulation of a Ross-Macdonald reflected SDE

# Load libraries

using OrdinaryDiffEq
using StochasticDiffEq
# using Plots
using DataFrames
# using DifferentialEquations.EnsembleAnalysis
using Main.Threads
using CSV
using ProgressBars
using Statistics

# using DiffEqMonteCarlo

# [] Later set this up as a function to parameters can be varied

# Set simulation parameters 

const deltat = 0.1f0 # timestep for Euler-Maruyama (EM) algorithm
const maxtime = 3650.0f0 # final time point
const time_vec = 1:deltat:maxtime # vector of timepoints at which to evaluate SDE
const num_timesteps = length(time_vec) # total number of timesteps for EM algorithm
const timespan = (0.0f0, maxtime)
const noise_rate_prototype = [0.0f0 0.0f0 0.0f0; 0.0f0 0.0f0 0.0f0] # prototype for noise

# Define the parameters
const b = 0.3f0 # biting rate
const Tvh = 0.5f0 # to vector transmission probability
const Nh = 10_000f0 # total number of humans
const Nv = 100_000f0 # total number of vectors
const muv = 0.1f0 # vector mortality rate
const gammah = 0.1f0 # human recovery rate
const alphab = 1 # toggle: environmental stochasticity in b
const alphat = 1 # toggle: environmental stochasticity in Tvh
const alpham = 1 # toggle: environmental stochasticity in muv
const q = [b, Tvh, Nh, Nv, muv, gammah]

## Initial conditions
const H0 = 0f0 # initial infected humans
const V0 = 10f0 # initial infected mosquitoes
const u0 = [H0, V0]

# Calculate Thv values from fixed R0 values
function Thv_from_R0(q, R0) 
    (b, Tvh, Nh, Nv, muv, gammah) = q
    Thv = (R0.^2) / ((b.^2 * Tvh * Nv) / (Nh * gammah * muv))
    return(Thv)
end

# Calculate endemic equilibrium values for deterministic case
function end_eqs(q, R0)
    # Set results to 0 where R0 < 1.0 using logical indexing
    if R0 .< 1.0f0
        V_end = 0.0f0
        H_end = 0.0f0
    else
    (b, Tvh, Nh, Nv, muv, gammah) = q
    Thv = Thv_from_R0(q, R0)
    rh = (Nv * b * Thv) / (Nh * gammah)
    rv = b * Tvh / muv
    
    H_end = Nh * (rh * rv - 1) / (rv + rh * rv)
    V_end = Nv * (rh * rv - 1) / (rh + rh * rv)

    # V_end = Nv * (R0 - 1) / (rh + R0) #Nv * (rh / (1 + rh)) * (R0 - 1) / R0
    # H_end = Nh * (rv / (1 + rv)) * (R0 - 1) / R0#Nh * (1 - (1 / (Nh^2 * rv))) / (1 + rh) 
    end
    end_vec = [H_end, V_end]
    return(end_vec)
end

# Initialize dataframes, if necessary


# Define the Ross-Macdonald Reflected SDE 

## Deterministic part
function dF_det!(du,u,p,t)
    # Name the variables
    (H,V) = u
    (Thv, sigma) = p #, b, Tvh, Nh, Nv, muv, gammah, alphab, alphat, alpham) = p

    du[1] = @views V * (b * Thv * (Nh - H) / Nh) - gammah * H
    du[2] = @views H * b * (Tvh / Nh) * (Nv - V) - muv * V
    return(du)
end


# without demographic stochasticity
function dF_det_no_demo!(du,u,p,t)
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
    du[1,1] = @views sqrt(max(0.0f0, (b * Thv * (Nh - H) / Nh) * V + gammah * H))
    du[2,1] = 0
    du[1,2] = 0
    du[2,2] = @views sqrt(max(0.0f0, b * (Tvh / Nh) * H * (Nv - V) + muv * V))
    # Environmental stochasticity
    du[1,3] = @views sigma * b * (Thv / Nh) * V * (Nh - H)
    du[2,3] = @views sigma * b * (Tvh / Nh) * H * (Nv - V) + sigma * b  * (Tvh / Nh) * H * (Nv - V) + sigma * muv *V
    return(du)
end

# without demographic stochasticity
function dF_stoch_no_demo!(du,u,p,t)
    # Name the variables
    (H,V) = u
    (Thv, sigma) = p #, b, Tvh, Nh, Nv, muv, gammah, alphab, alphat, alpham) = p
    # Demographic stochasticity
    du[1,1] = 0f0
    du[2,1] = 0f0
    du[1,2] = 0f0
    du[2,2] = 0f0
    # Environmental stochasticity
    du[1,3] = @views sigma * b * (Thv / Nh) * V * (Nh - H)
    du[2,3] = @views sigma * b * (Tvh / Nh) * H * (Nv - V) + sigma * b  * (Tvh / Nh) * H * (Nv - V) + sigma * muv *V
    return(du)
end

## Reflections
# --- In Julia these are done via callbacks
condition_H(u,t,integrator) = true
function affect_H!(integrator)
    # Reflect if H goes below zero
    if integrator.u[1] < 0.0f0
        integrator.u[1] = 0.0f0
        # Reflect if H goes above Nh
    elseif integrator.u[1] > Nh
        integrator.u[1] = Nh
    end
end
cb_H = DiscreteCallback(condition_H, affect_H!; save_positions = (false, false))

condition_V(u,t,integrator) = true
function affect_V!(integrator)
    # Reflect if V goes below zero
    if integrator.u[2] < 0.0f0
        integrator.u[2] = 0.0f0
        # Reflect if V goes above Nv
    elseif integrator.u[2] > Nv
        integrator.u[2] = Nv
    end
end
cb_V = DiscreteCallback(condition_V, affect_V!; save_positions = (false, false))

### Stop simulations if case counts are low enough
function condition_terminate(u, t, integrator)
    u[2] < 1.0f0 && u[1] < 1.0f0
end
affect_terminate!(integrator) = terminate!(integrator)
cb_terminate = DiscreteCallback(condition_terminate, affect_terminate!; save_positions = (false, false))

## Collect callbacks
# cbs = CallbackSet(cb_H0, cb_HNh, cb_V0, cb_VNv, cb_terminate)
cbs = CallbackSet(cb_H, cb_V, cb_terminate)

# Define parameter values to iterate over
const R0s = 0f0:0.05f0:5f0 #0f0:0.025f0:5f0
Thvs = Thv_from_R0(q, R0s) # used to vary R0
const sigmas = 0f0:0.05f0:2f0 # levels of environmental noise

parameter_values = [(Thv, sigma) for Thv in Thvs, sigma in sigmas]


# Get example trajectories for grid plot ----------------------------------------------------------------------------
# Initialize an empty DataFrame to store results
param_names = ["H", "V"]

function run_sims(det_equations, stoch_equations, num_runs, parameter_values)
    results = DataFrame(time = Float64[], Thv = Float64[], sigma = Float64[], run = Int[], H = Float64[], V = Float64[])

    for i in ProgressBar(eachindex(parameter_values))

        ## Set up ensemble SDE problem
        prob = SDEProblem(det_equations, stoch_equations, u0, timespan, parameter_values[i], noise_rate_prototype = noise_rate_prototype, callback = cbs)
        ensembleprob = EnsembleProblem(prob)
        ## Run SDE solver
        sol = solve(ensembleprob, EM(), dt = 0.1f0, EnsembleThreads(); trajectories = num_runs, saveat = 1.0f0)

        # Collect results into a tidy DataFrame
        for run_id in 1:num_runs
            trajectory = sol[run_id]  # Get the individual trajectory solution
            # Because the callbacks save at each time step before reflection, we have to remove the points 
            # that are biologically unreasonable first
            trajectory = unique([(trajectory.t[i], trajectory.u[i][1], trajectory.u[i][2]) 
                    for i in eachindex(trajectory.t) 
                    if 0 ≤ trajectory.u[i][1] ≤ Nh && 0 ≤ trajectory.u[i][2] ≤ Nv])
            times = first.(trajectory)

            results_append = DataFrame(
                time = first.(trajectory),
                Thv = fill(parameter_values[i][1], length(times)),
                sigma = fill(parameter_values[i][2], length(times)),
                run = fill(run_id, length(times)),
                H = getindex.(trajectory,2),
                V = last.(trajectory)
            )
            
            # Append to the main results DataFrame
            append!(results, results_append)
        end
    end
    return results
end

# Collect summarized outputs ----------------------------------------------------------------------------
prob = SDEProblem(dF_det!, dF_stoch!, u0, timespan, parameter_values[100], noise_rate_prototype = noise_rate_prototype, callback = cbs)
ensembleprob = EnsembleProblem(prob)
## Run SDE solver

# Precompile the EM solver
sol = solve(ensembleprob, EM(), dt = 0.1f0, EnsembleThreads(); trajectories = 10::Int)

# Collect outputs from SDE simulations
function collect_outputs(det_equations, stoch_equations, num_runs, parameter_values)
    # Preallocate a DataFrame to store results for each trajectory and parameter combination
    results_init = DataFrame(
        Thv = Float32[], sigma = Float32[], name = String[], statistic = String[], value = Float32[])
    results = results_init

    for i in ProgressBar(eachindex(parameter_values))
        # Get trajectories
        prob = SDEProblem(det_equations, stoch_equations, u0, timespan, parameter_values[i], noise_rate_prototype = noise_rate_prototype, callback = cbs)
        ensembleprob = EnsembleProblem(prob)
        ## Run SDE solver
        sol = solve(ensembleprob, EM(), dt = 0.1f0, EnsembleSplitThreads(); trajectories = num_runs, saveat = 1.0f0)

        # Thread-local storage for results
        thread_results = Vector{NamedTuple{(:Thv, :sigma, :run, :max_value, :duration,:peak_time, :exceeded_10, :exceeded_100, :positive_at_final, :positive_duration, :zero_cases,:duration_dieout),Tuple{Float32, Float32, Int, Float32, Float32, Float32, Bool, Bool, Bool, Float32, Bool, Union{Missing,Float32}}}}(undef, num_runs)

        Thv, sigma = parameter_values[i]
        
        

        # @threads for run_id in ProgressBar(1:num_runs)
        #     trajectory = sol[run_id]  # Get the individual trajectory solution
        #     nt = length(trajectory.t)
        
        #     # Pre-allocate H_values
        #     H_values = Vector{Float32}(undef, nt)
        #     @inbounds for i in 1:nt
        #         H_values[i] = trajectory.u[i][1]
        #     end
        
        #     max_value = maximum(H_values)
        
        #     # Compute valid_times using an explicit loop (avoid generator overhead)
        #     valid_times = Float32[]
        #     for i in 1:nt
        #         @inbounds begin
        #             if trajectory.u[i][1] > 1.0f0 && trajectory.u[i][2] > 1.0f0
        #                 push!(valid_times, trajectory.t[i])
        #             end
        #         end
        #     end
        #     duration = isempty(valid_times) ? 0.0f0 : maximum(valid_times)/365.0f0
        
        #     # peak_time is taken from trajectory.t at the index where the second state is maximum.
        #     # (Adjust if you intended something else.)
        #     idx = argmax(trajectory[2, :])  # Assumes trajectory[2, :] returns the second state variable over time.
        #     peak_time = (max_value < 1.0f0) ? 0.0f0 : trajectory.t[idx]
        
        #     exceeded_10 = false
        #     exceeded_100 = false
        #     @inbounds for v in H_values
        #         if v > 10.0f0
        #             exceeded_10 = true
        #         end
        #         if v > 100.0f0
        #             exceeded_100 = true
        #         end
        #     end
        
        #     positive_at_final = duration > 9.999f0 && H_values[end] > 1.0f0
        
        #     # Compute positive_duration as the total time duration where H > 1 in consecutive time intervals
        #     positive_duration = 0.0f0
        #     @inbounds for i in 1:nt-1
        #         if trajectory.u[i][1] > 1.0f0 && trajectory.u[i+1][1] > 1.0f0
        #             positive_duration += trajectory.t[i+1] - trajectory.t[i]
        #         end
        #     end
        #     positive_duration /= 365.0f0
        
        #     zero_cases = false
        #     @inbounds for i in 2:nt
        #         if H_values[i] < 1.0f0
        #             zero_cases = true
        #             break
        #         end
        #     end
        
        #     duration_dieout = positive_at_final ? missing : duration
        
        #     # Store the results in the preallocated array.
        #     thread_results[run_id] = (Thv, sigma, run_id, max_value, duration, peak_time, 
        #                               exceeded_10, exceeded_100, positive_at_final, positive_duration,
        #                               zero_cases, duration_dieout)
        # end

        # # Analyze each trajectory
        # @threads for run_id in ProgressBar(1:num_runs)
        #     trajectory = sol[run_id]  # Get the individual trajectory solution

        #     # Extract the values and times for the state variables
        #     H_values = [u[1] for u in @views trajectory.u]

        #     # Calculate the required statistics
        #     # Maximum case count
        #     max_value = maximum(@views H_values)
        #     # Get times where there is at least one infection in each population
        #     valid_times = (trajectory.t[i] for i in eachindex(trajectory.t) if trajectory.u[i][1] > 1.0f0 && trajectory.u[i][2] > 1.0f0)
        #     # Calculate the latest date where there are cases in both hosts and vectors
        #     duration = isempty(valid_times) ? 0.0f0 : maximum(valid_times)/365.0f0
        #     peak_time = ifelse(max_value < 1.0f0, 0.0f0, trajectory.t[argmax(trajectory[2, :])])
        #     # Check if host case counts ever exceeded 10 or 100
        #     exceeded_10 = any(v > 10.0f0 for v in H_values)
        #     exceeded_100 = any(v > 100.0f0 for v in H_values)
        #     # Check whether infections lasted all ten years
        #     positive_at_final = duration > 9.999f0 && @views H_values[end] > 1.0f0
        #     # Calculate positive_duration (outbreak length)
        #     valid_pairs = [(trajectory.t[i], trajectory.u[i][1], trajectory.t[i+1], trajectory.u[i+1][1]) for i in 1:length(trajectory.t)-1 if trajectory.u[i][1] > 1.0f0 && trajectory.u[i+1][1] > 1.0f0]
        #     positive_duration = isempty(valid_pairs) ? 0.0f0 : sum(t2 - t1 for (t1, u1, t2, u1_next) in valid_pairs)/365.0f0
        #     # See if cases ever dropped to zero in hosts (besides at the initial time point)
        #     zero_cases = any(@views H_values[2:end] .< 1.0f0)
        #     duration_dieout = positive_at_final == true ? missing : duration

        #     # Store results in thread-local array
        #     thread_results[run_id] = (Thv, sigma, run_id, max_value, duration, peak_time, exceeded_10, exceeded_100, positive_at_final, positive_duration, zero_cases, duration_dieout)
        # end

        @threads for run_id in 1:num_runs
            traj = sol[run_id]              # Get the trajectory (assumed independent per run)
            nt = length(traj.t)
            @inbounds begin
                # Initialize local accumulators; use Float32 for consistency.
                max_val       = -typemax(Float32)
                max_H_i       = -typemax(Float32)
                exp10         = false
                exp100        = false
                valid_time    = 0.0f0         # We will take maximum time where both populations > 1.
                pos_duration  = 0.0f0         # Sum of positive intervals (in same time units as traj.t)
                zero_cases    = false
        
                # Loop over trajectory indices once.
                for i in 1:nt
                    # Extract the host (H) count and second state from the state vector.
                    H_i    = traj.u[i][1]
                    V_i  = traj.u[i][2]
                    t_i    = traj.t[i]
                    
                    # Update max H.
                    if H_i > max_val
                        max_val = H_i
                        max_H_i = i
                    end
        
                    # Check for exceeded thresholds.
                    if H_i > 10.0f0
                        exp10 = true
                    end
                    if H_i > 100.0f0
                        exp100 = true
                    end
        
                    # Update valid_time if both H and V exceed 1.
                    if (traj.u[i][1] > 1.0f0) && (traj.u[i][2] > 1.0f0) && (t_i > valid_time)
                        valid_time = t_i
                    end
                    # Accumulate positive duration for intervals where H is > 1.
                    if i < nt
                        if traj.u[i][1] > 1.0f0 && traj.u[i+1][1] > 1.0f0
                            pos_duration += 1
                        end
                    end
        
                    # Check for zero cases (after the initial time point)
                    if i > 1 && H_i < 1.0f0
                        zero_cases = true
                    end
                end  # end of trajectory loop
        
                # Calculate duration, peak time, etc.
                duration          = valid_time == 0.0f0 ? 0.0f0 : valid_time / 365.0f0
                peak_time         = (max_val < 1.0f0) ? 0.0f0 : traj.t[max_H_i]
                positive_at_final = (duration > 9.999f0) && (traj.u[nt][1] > 1.0f0)
                pos_duration      = pos_duration / 365.0f0
                duration_dieout   = positive_at_final ? missing : duration
        
                # Store results into the preallocated array.
                thread_results[run_id] = (Thv, sigma, run_id, max_val, duration, peak_time,
                                          exp10, exp100, positive_at_final, pos_duration,
                                          zero_cases, duration_dieout)
            end
        end

        thread_results = DataFrame(thread_results)
    
        # Thread-local storage (prevents concurrent modification)
        temp_results = DataFrame(
            Thv = Float32[], sigma = Float32[], name = String[], mean = Float32[], median = Float32[], variance = Float32[], q_25 = Float32[], q_75 = Float32[])

            for col in [:max_value, :duration, :peak_time, :exceeded_10, :exceeded_100, :positive_at_final, :positive_duration, :zero_cases, :duration_dieout]
            col_values = skipmissing(thread_results[!, col])  # Extract column values
            if isempty(col_values)
                push!(temp_results, (Thv, sigma, string(col), NaN32, NaN32, NaN32, NaN32, NaN32, ))
            else
                push!(temp_results, (Thv, sigma, string(col), mean(col_values), median(col_values), var(col_values), quantile(col_values, 0.25), quantile(col_values, 0.75)))
            end
        end
        temp_results_long = stack(DataFrame(temp_results), [:mean, :median, :variance, :q_25, :q_75], variable_name=:statistic, value_name=:value)

        # Append results safely after threading
        append!(results, temp_results_long)
    end
    return results
end


# Collect outputs from SDE simulations
function raw_outputs(det_equations, stoch_equations, num_runs, parameter_values)
    # Preallocate a DataFrame to store results for each trajectory and parameter combination
    # results_init = DataFrame(
        # Thv = Float32[], sigma = Float32[], run = Float32[], statistic = String[], value = Float32[])
    results_init = Vector{NamedTuple{(:Thv, :sigma, :run, :max_value, :max_time,:exceeded_10, :exceeded_100, :positive_at_final, :positive_duration, :zero_cases,:duration_dieout),Tuple{Float32, Float32, Int, Float32, Float32, Bool, Bool, Bool, Float32, Bool, Float32}}}(undef, 0)
    results = DataFrame(results_init)

    for i in ProgressBar(eachindex(parameter_values))
        # Get trajectories
        prob = SDEProblem(det_equations, stoch_equations, u0, timespan, parameter_values[i], noise_rate_prototype = noise_rate_prototype, callback = cbs)
        ensembleprob = EnsembleProblem(prob)
        ## Run SDE solver
        sol = solve(ensembleprob, EM(), dt = 0.1f0, EnsembleSplitThreads(); trajectories = num_runs, saveat = 1.0f0, save_everystep = false)
        # Thread-local storage for results
        thread_results = Vector{NamedTuple{(:Thv, :sigma, :run, :max_value, :max_time,:exceeded_10, :exceeded_100, :positive_at_final, :positive_duration, :zero_cases,:duration_dieout),Tuple{Float32, Float32, Int, Float32, Float32, Bool, Bool, Bool, Float32, Bool, Float32}}}(undef, num_runs)

        Thv = parameter_values[i][1]
        sigma = parameter_values[i][2]
        # Analyze each trajectory
        @threads for run_id in 1:num_runs
            trajectory = sol[run_id]  # Get the individual trajectory solution

            # Because the callbacks save at each time step before reflection, we have to remove the points 
            # that are biologically unreasonable first
            filtered_trajectory = unique([(t, u1, u2) for (t, (u1, u2)) in zip(trajectory.t, trajectory.u) if 0 <= u1 <= Nh && 0 <= u2 <= Nv])
            
            # Extract the values and times for the state variables
            H_values = getindex.(filtered_trajectory, 2)

            # Calculate the required statistics
            # Maximum case count
            max_value = maximum(H_values) # only needs largest H value
            # Get times where there is at least one infection in each population
            valid_times = (trajectory.t[i] for i in eachindex(trajectory.t) if trajectory.u[i][1] > 1.0f0 && trajectory.u[i][2] > 1.0f0)
            # Calculate the latest date where there are cases in both hosts and vectors
            max_time = isempty(valid_times) ? 0.0f0 : maximum(valid_times)/365.0f0
            # Check if host case counts ever exceeded 10 or 100
            exceeded_10 = any(v > 10.0f0 for v in H_values)
            exceeded_100 = any(v > 100.0f0 for v in H_values)
            # Check whether infections lasted all ten years
            positive_at_final = max_time > 9.99f0 && H_values[end] > 1.0f0
            # Calculate positive_duration (outbreak length)
            valid_pairs = [(t1, u1, t2, u1_next) for ((t1, u1, _), (t2, u1_next, _)) in zip(filtered_trajectory[1:end-1], filtered_trajectory[2:end]) if u1 > 1.0f0 && u1_next > 1.0f0]
            positive_duration = isempty(valid_pairs) ? 0.0f0 : sum(t2 - t1 for (t1, u1, t2, u1_next) in valid_pairs)/365.0f0
            # See if cases ever dropped to zero in hosts (besides at the initial time point)
            zero_cases = any(H_values[2:end] .< 1.0f0)
            duration_dieout = positive_at_final == true ? NaN : max_time
            # [] add dieout duration

            # Store results in thread-local array
            thread_results[run_id] = (Thv, sigma, run_id, max_value, max_time, exceeded_10, exceeded_100, positive_at_final, positive_duration, zero_cases, duration_dieout)
            # Append results for this trajectory to the DataFrame
            # push!(results, (Thv, sigma, run_id, max_value, max_time, exceeded_10, exceeded_100, positive_at_final, positive_duration))
        end
        thread_results = DataFrame(thread_results)

        # Append results safely after threading
        append!(results, thread_results)
    end
    return results
end

# Collect outputs from ODE simulations
function collect_outputs_det(det_equations, parameter_values)
    
    # Initialize a DataFrame to store results for each trajectory and parameter combination
    results = DataFrame(Thv = Float32[], R0 = Float32[], max_value = Float32[], max_time = Float32[], exceeded_10 = Bool[],
    exceeded_100 = Bool[], positive_at_final = Bool[], positive_duration = Float32[], H_end_exact = Float32[], H_end_num = Float32[], 
    V_end_exact = Float32[], V_end_num = Float32[])

    for i in ProgressBar(eachindex(parameter_values))
        R0 = parameter_values[i]
        Thv = Thv_from_R0(q, R0)
        # Get trajectories
        prob = ODEProblem(det_equations, u0, timespan, [Thv, 0.0f0]) # timespan
        ## Run SDE solver
        sol = solve(prob; saveat = 1.0f0)

        # Analyze each trajectory
        H_end_exact, V_end_exact = end_eqs(q, R0)

        trajectory = sol  # Get the individual trajectory solution

        # Because the callbacks save at each time step before reflection, we have to remove the points 
        # that are biologically unreasonable first
        filtered_trajectory = unique([(t, u1, u2) for (t, (u1, u2)) in zip(trajectory.t, trajectory.u) if 0 <= u1 <= Nh && 0 <= u2 <= Nv])
        
        # Extract the values and times for the state variables
        H_values = getindex.(filtered_trajectory, 2)
        V_values = getindex.(filtered_trajectory, 3)
        H_end_num = H_values[end]
        V_end_num = V_values[end]

        # Calculate the required statistics
        # Maximum case count
        max_value = maximum(H_values)
        # Get times where there is at least one infection in each population
        valid_times = (trajectory.t[i] for i in eachindex(trajectory.t) if trajectory.u[i][1] > 1.0f0 && trajectory.u[i][2] > 1.0f0)
        # Calculate the latest date where there are cases in both hosts and vectors
        max_time = isempty(valid_times) ? 0.0f0 : maximum(valid_times)
        # Check if host case counts ever exceeded 10 or 100
        exceeded_10 = any(v > 10.0f0 for v in H_values)
        exceeded_100 = any(v > 100.0f0 for v in H_values)
        # Check whether infections lasted all ten years
        positive_at_final = max_time == maximum(trajectory.t) && H_values[end] > 1.0f0
        # Calculate positive_duration (outbreak length)
        valid_pairs = [(t1, u1, t2, u1_next) for ((t1, u1, _), (t2, u1_next, _)) in zip(filtered_trajectory[1:end-1], filtered_trajectory[2:end]) if u1 > 1.0f0 && u1_next > 1.0f0]
        positive_duration = isempty(valid_pairs) ? 0.0f0 : sum(t2 - t1 for (t1, u1, t2, u1_next) in valid_pairs)
        
        Thv = parameter_values[1]
        

        # Append results for this trajectory to the DataFrame
        push!(results, (Thv, parameter_values[i], max_value, max_time, exceeded_10, exceeded_100, positive_at_final, positive_duration, H_end_exact, H_end_num, V_end_exact, V_end_num))

    end
    return(results)
end
