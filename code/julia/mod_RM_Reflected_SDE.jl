# Numerical simulation of a Ross-Macdonald reflected SDE

# Load libraries
using OrdinaryDiffEq
using StochasticDiffEq
using DataFrames
using Main.Threads
using CSV
using ProgressBars
using Statistics

# using DiffEqMonteCarlo

# [] Later set this up as a function so parameters can be varied

# Set simulation parameters 

const deltat = 0.1f0 # timestep for Euler-Maruyama (EM) algorithm
const maxtime = 3650.0f0 # final time point
const time_vec = 1:deltat:maxtime # vector of timepoints at which to evaluate SDE
const num_timesteps = length(time_vec) # total number of timesteps for EM algorithm
const timespan = (0.0f0, maxtime)
const noise_rate_prototype = [0.0f0 0.0f0 0.0f0; 0.0f0 0.0f0 0.0f0] # prototype for noise

# Define the constant parameters
const muv = 0.1f0 # vector mortality rate
const gammah = 0.1f0 # human recovery rate
const alphab = 1 # toggle: environmental stochasticity in b
const alphat = 1 # toggle: environmental stochasticity in Tvh
const alpham = 1 # toggle: environmental stochasticity in muv
# Note: b, Tvh, Nh, Nv are now variable parameters passed in the parameter tuple

## Initial conditions
const H0 = 0f0 # initial infected humans
const V0 = 10f0 # initial infected mosquitoes
const u0 = [H0, V0]

# Calculate Thv values from fixed R0 values
function Thv_from_R0(b, Tvh, Nh, Nv, R0) 
    Thv = (R0.^2) / ((b.^2 * Tvh * Nv) / (Nh * gammah * muv))
    return(Thv)
end

# Calculate R0 values from parameters
function R0_calc(b, Thv, Tvh, Nh, Nv) 
    R0 = sqrt((b.^2 * Thv * Tvh * Nv) / (Nh * gammah * muv))
    return(R0)
end

# Calculate endemic equilibrium values for deterministic case
function end_eqs(b, Tvh, Nh, Nv, R0)
    # Set results to 0 where R0 < 1.0 using logical indexing
    if R0 .< 1.0f0
        V_end = 0.0f0
        H_end = 0.0f0
    else
    Thv = Thv_from_R0(b, Tvh, Nh, Nv, R0)
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

# Define the Ross-Macdonald Reflected SDE 

## Deterministic part
function dF_det!(du,u,p,t)
    # Name the variables
    (H,V) = u
    (b, Tvh, Nh, Nv, Thv, sigma) = p

    du[1] = @views V * (b * Thv * (Nh - H) / Nh) - gammah * H
    du[2] = @views H * b * (Tvh / Nh) * (Nv - V) - muv * V
    return(du)
end


# without demographic stochasticity
function dF_det_no_demo!(du,u,p,t)
    # Name the variables
    (H,V) = u
    (b, Tvh, Nh, Nv, Thv, sigma) = p

    du[1] = @views (b * Thv * (Nh - H) / Nh) * V - gammah * H
    du[2] = @views b * (Tvh / Nh) * H * (Nv - V) - muv * V
    return(du)
end

## Stochastic part
function dF_stoch!(du,u,p,t)
    # Name the variables
    (H,V) = u
    (b, Tvh, Nh, Nv, Thv, sigma) = p
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
    (b, Tvh, Nh, Nv, Thv, sigma) = p
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
    elseif integrator.u[1] > integrator.p[3]  # Nh is the 3rd parameter
        integrator.u[1] = integrator.p[3]
    end
end
cb_H = DiscreteCallback(condition_H, affect_H!; save_positions = (false, false))

condition_V(u,t,integrator) = true
function affect_V!(integrator)
    # Reflect if V goes below zero
    if integrator.u[2] < 0.0f0
        integrator.u[2] = 0.0f0
        # Reflect if V goes above Nv
    elseif integrator.u[2] > integrator.p[4]  # Nv is the 4th parameter
        integrator.u[2] = integrator.p[4]
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

# Define default values for the variable parameters
const b_default = 0.3f0
const Tvh_default = 0.5f0
const Nh_default = 10_000f0
const Nv_default = 100_000f0

# Define parameter values to iterate over
const R0s = 0f0:0.05f0:5f0 #0f0:0.025f0:5f0
const sigmas = 0f0:0.05f0:2f0 # levels of environmental noise

# Create parameter tuple as (b, Tvh, Nh, Nv, Thv, sigma)
parameter_values = [(b_default, Tvh_default, Nh_default, Nv_default, Thv_from_R0(b_default, Tvh_default, Nh_default, Nv_default, R0), sigma) 
                    for R0 in R0s, sigma in sigmas]

# Get example trajectories for grid plot ----------------------------------------------------------------------------
# Initialize an empty DataFrame to store results
param_names = ["H", "V"]

function run_sims(det_equations, stoch_equations, num_runs, parameter_values)
    results = DataFrame(time = Float64[], b = Float32[], Tvh = Float32[], Nh = Float32[], Nv = Float32[], Thv = Float32[], sigma = Float32[], run = Int[], H = Float64[], V = Float64[])

    for i in ProgressBar(eachindex(parameter_values))
        b, Tvh, Nh, Nv, Thv, sigma = parameter_values[i]

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
                b = fill(b, length(times)),
                Tvh = fill(Tvh, length(times)),
                Nh = fill(Nh, length(times)),
                Nv = fill(Nv, length(times)),
                Thv = fill(Thv, length(times)),
                sigma = fill(sigma, length(times)),
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

# Collect outputs from SDE simulations (optimized for 100K+ runs)
function collect_outputs(det_equations, stoch_equations, num_runs, parameter_values; output_file="collect_outputs.csv")
    # Pre-open file for writing
    col_names = [:max_value, :duration, :peak_time, :exceeded_10, :exceeded_100, :positive_at_final, :positive_duration, :zero_cases, :duration_dieout]
    first_write = true
    
    for i in ProgressBar(eachindex(parameter_values))
        # Get trajectories
        prob = SDEProblem(det_equations, stoch_equations, u0, timespan, parameter_values[i], noise_rate_prototype = noise_rate_prototype, callback = cbs)
        ensembleprob = EnsembleProblem(prob)
        ## Run SDE solver
        sol = solve(ensembleprob, EM(), dt = 0.1f0, EnsembleSplitThreads(); trajectories = num_runs, saveat = 1.0f0, verbose = false, progress = false)

        # Thread-local storage for results
        thread_results = Vector{NamedTuple{(:b, :Tvh, :Nh, :Nv, :Thv, :sigma, :run, :max_value, :duration,:peak_time, :exceeded_10, :exceeded_100, :positive_at_final, :positive_duration, :zero_cases,:duration_dieout),Tuple{Float32, Float32, Float32, Float32, Float32, Float32, Int, Float32, Float32, Float32, Bool, Bool, Bool, Float32, Bool, Union{Missing,Float32}}}}(undef, num_runs)

        b, Tvh, Nh, Nv, Thv, sigma = parameter_values[i]

        @threads for run_id in 1:num_runs
            traj = sol[run_id]
            nt = length(traj.t)
            @inbounds begin
                max_val = -typemax(Float32)
                max_H_i = -typemax(Float32)
                exp10 = false
                exp100 = false
                valid_time = 0.0f0
                pos_duration = 0.0f0
                zero_cases = false
        
                for j in 1:nt
                    H_j = traj.u[j][1]
                    V_j = traj.u[j][2]
                    t_j = traj.t[j]
                    
                    if H_j > max_val
                        max_val = H_j
                        max_H_i = j
                    end
                    if H_j > 10.0f0
                        exp10 = true
                    end
                    if H_j > 100.0f0
                        exp100 = true
                    end
                    if (H_j > 1.0f0) && (V_j > 1.0f0) && (t_j > valid_time)
                        valid_time = t_j
                    end
                    if j < nt && traj.u[j][1] > 1.0f0 && traj.u[j+1][1] > 1.0f0
                        pos_duration += 1
                    end
                    if j > 1 && H_j < 1.0f0
                        zero_cases = true
                    end
                end
        
                duration = valid_time == 0.0f0 ? 0.0f0 : valid_time / 365.0f0
                peak_time = (max_val < 1.0f0) ? 0.0f0 : traj.t[max_H_i]
                positive_at_final = (duration > 9.999f0) && (traj.u[nt][1] > 1.0f0)
                pos_duration = pos_duration / 365.0f0
                duration_dieout = positive_at_final ? missing : duration
        
                thread_results[run_id] = (b, Tvh, Nh, Nv, Thv, sigma, run_id, max_val, duration, peak_time, exp10, exp100, positive_at_final, pos_duration, zero_cases, duration_dieout)
            end
        end

        # Convert to DataFrame and compute stats in one shot
        df_raw = DataFrame(thread_results)
        batch_results = DataFrame(b = Float32[], Tvh = Float32[], Nh = Float32[], Nv = Float32[], Thv = Float32[], sigma = Float32[], 
                                 name = String[], statistic = String[], value = Float32[])
        
        # Compute all statistics in a single pass per column
        for col_name in col_names
            col_values = skipmissing(df_raw[!, col_name])
            
            if isempty(col_values)
                for stat in [:mean, :median, :variance, :q_25, :q_75, :min, :max]
                    push!(batch_results, (b, Tvh, Nh, Nv, Thv, sigma, string(col_name), string(stat), NaN32))
                end
            else
                # Convert Bool to Float32 if needed
                if eltype(col_values) == Bool
                    col_values = Float32.(col_values)
                end
                
                m = mean(col_values)
                med = median(col_values)
                v = var(col_values)
                q25 = quantile(col_values, 0.25f0)
                q75 = quantile(col_values, 0.75f0)
                min = minimum(col_values)
                max = maximum(col_values)
                
                push!(batch_results, (b, Tvh, Nh, Nv, Thv, sigma, string(col_name), "mean", Float32(m)))
                push!(batch_results, (b, Tvh, Nh, Nv, Thv, sigma, string(col_name), "median", Float32(med)))
                push!(batch_results, (b, Tvh, Nh, Nv, Thv, sigma, string(col_name), "variance", Float32(v)))
                push!(batch_results, (b, Tvh, Nh, Nv, Thv, sigma, string(col_name), "q_25", Float32(q25)))
                push!(batch_results, (b, Tvh, Nh, Nv, Thv, sigma, string(col_name), "q_75", Float32(q75)))
                push!(batch_results, (b, Tvh, Nh, Nv, Thv, sigma, string(col_name), "min", Float32(min)))
                push!(batch_results, (b, Tvh, Nh, Nv, Thv, sigma, string(col_name), "max", Float32(max)))
            end
        end
        
        # Write batch to CSV
        CSV.write(output_file, batch_results, append=!first_write)
        first_write = false
    end
    
    # Return filename
    return output_file
end


# Collect outputs from SDE simulations (optimized for large parameter spaces)
function raw_outputs(det_equations, stoch_equations, num_runs, parameter_values; output_file="raw_outputs.csv")
    # Write to CSV incrementally instead of building giant DataFrame
    first_write = true
    
    for i in ProgressBar(eachindex(parameter_values))
        # Get trajectories
        prob = SDEProblem(det_equations, stoch_equations, u0, timespan, parameter_values[i], noise_rate_prototype = noise_rate_prototype, callback = cbs)
        ensembleprob = EnsembleProblem(prob)
        ## Run SDE solver
        sol = solve(ensembleprob, EM(), dt = 0.1f0, EnsembleSplitThreads(); trajectories = num_runs, saveat = 1.0f0, save_everystep = false)
        
        # Pre-allocate batch results
        b, Tvh, Nh, Nv, Thv, sigma = parameter_values[i]
        batch_results = Vector{NamedTuple{(:b, :Tvh, :Nh, :Nv, :Thv, :sigma, :run, :max_value, :max_time,:exceeded_10, :exceeded_100, :positive_at_final, :positive_duration, :zero_cases,:duration_dieout),Tuple{Float32, Float32, Float32, Float32, Float32, Float32, Int, Float32, Float32, Bool, Bool, Bool, Float32, Bool, Float32}}}(undef, num_runs)

        # Analyze each trajectory with single-pass optimization
        @threads for run_id in 1:num_runs
            trajectory = sol[run_id]
            nt = length(trajectory.t)
            
            @inbounds begin
                # Single pass through trajectory computing all statistics
                max_value = 0.0f0
                last_valid_time = 0.0f0
                exceeded_10 = false
                exceeded_100 = false
                pos_duration = 0.0f0
                zero_cases = false
                prev_H_valid = false
                
                for j in 1:nt
                    H_j = trajectory.u[j][1]
                    V_j = trajectory.u[j][2]
                    t_j = trajectory.t[j]
                    
                    # Skip biologically unreasonable points
                    if H_j < 0.0f0 || H_j > Nh || V_j < 0.0f0 || V_j > Nv
                        prev_H_valid = false
                        continue
                    end
                    
                    # Update max value
                    if H_j > max_value
                        max_value = H_j
                    end
                    
                    # Check thresholds
                    if H_j > 10.0f0
                        exceeded_10 = true
                    end
                    if H_j > 100.0f0
                        exceeded_100 = true
                    end
                    
                    # Track valid time (both H and V > 1)
                    if H_j > 1.0f0 && V_j > 1.0f0
                        last_valid_time = t_j
                    end
                    
                    # Accumulate positive duration
                    if H_j > 1.0f0
                        if prev_H_valid
                            pos_duration += 1.0f0  # Assume 1 day timestep, will divide by 365 later
                        end
                        prev_H_valid = true
                    else
                        if j > 1
                            zero_cases = true
                        end
                        prev_H_valid = false
                    end
                end
                
                # Calculate final statistics
                max_time = last_valid_time == 0.0f0 ? 0.0f0 : last_valid_time / 365.0f0
                pos_duration = pos_duration / 365.0f0
                positive_at_final = (max_time > 9.99f0) && (trajectory.u[nt][1] > 1.0f0)
                duration_dieout = positive_at_final ? NaN32 : max_time
                
                # Store result
                batch_results[run_id] = (b, Tvh, Nh, Nv, Thv, sigma, Int32(run_id), max_value, max_time, exceeded_10, exceeded_100, positive_at_final, pos_duration, zero_cases, duration_dieout)
            end
        end
        
        # Convert to DataFrame and write
        batch_df = DataFrame(batch_results)
        CSV.write(output_file, batch_df, append=!first_write)
        first_write = false
    end
    
    # Return the full dataset
    return CSV.read(output_file, DataFrame)
end

# Collect outputs from ODE simulations
function collect_outputs_det(det_equations, R0_values)
    
    # Initialize a DataFrame to store results for each trajectory and parameter combination
    results = DataFrame(b = Float32[], Tvh = Float32[], Nh = Float32[], Nv = Float32[], R0 = Float32[], Thv = Float32[], max_value = Float32[], max_time = Float32[], exceeded_10 = Bool[],
    exceeded_100 = Bool[], positive_at_final = Bool[], positive_duration = Float32[], H_end_exact = Float32[], H_end_num = Float32[], 
    V_end_exact = Float32[], V_end_num = Float32[])

    for R0 in ProgressBar(R0_values)
        b, Tvh, Nh, Nv = b_default, Tvh_default, Nh_default, Nv_default
        Thv = Thv_from_R0(b, Tvh, Nh, Nv, R0)
        # Get trajectories
        prob = ODEProblem(det_equations, u0, timespan, [b, Tvh, Nh, Nv, Thv, 0.0f0]) # timespan
        ## Run ODE solver
        sol = solve(prob; saveat = 1.0f0)

        # Analyze each trajectory
        H_end_exact, V_end_exact = end_eqs(b, Tvh, Nh, Nv, R0)

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
        

        # Append results for this trajectory to the DataFrame
        push!(results, (b, Tvh, Nh, Nv, R0, Thv, max_value, max_time, exceeded_10, exceeded_100, positive_at_final, positive_duration, H_end_exact, H_end_num, V_end_exact, V_end_num))

    end
    return(results)
end
