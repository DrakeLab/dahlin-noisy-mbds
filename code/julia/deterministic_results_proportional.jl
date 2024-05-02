
# Load necessary packages_______________________________________________________
using CSV
using DataFrames
using DifferentialEquations
using Plots
using Statistics
# using PlotlyJS

# Define functions_____________________________________________________________
# base function for mosquito borne disease
function f(du,u,p,t)
  #(H,V) = u
  #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
  du[1] = @views (p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2] - p[6] * u[1]
  #du[1] = (b * τₕᵥ * (Nₕ - H) * V / Nₕ) - γₕ * H
  du[2] = @views (p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] - p[7] * u[2]
  #du[2] = (b * τᵥₕ * (Nᵥ - V) * H / Nₕ) - μᵥ * V]

  return(du)
end

# Function: Calculate tVH values from specified R0 values
function tHV_from_R0(p, R0) 
  (b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
  tHV = (R0.^2) / ((b.^2 * τᵥₕ * Nᵥ) / (Nₕ * γₕ * μᵥ))
  return(tHV)
end

# Define parameters___________________________________________________________________________
#Parameters: biting rate, THV, TVH, NH, NV, recovery, mortality, sigma, alphab, alpham, alphat
p = [0.3, 0.02450000, 0.5, 10000, 100000, 0.1, 0.1, 0.0125, 1.0, 1.0, 1.0]

#Rate of transmission from hosts to vectors-We varied this parameter to get different R0 values
#Code for calculating Tvh for a given R0 value is in the main R code


#initial conditions->10 infected mosquitos to start
u0 = [0.0, 10.0]
#10 year simulation with time step of 1 day
tspan = [0.0, 11*365]
dt = 1


# Get deterministic ODE solution trajectories______________________________________________
R0s = [0.95, 1.05, 1.25, 2, 4, 6.5]
THVs = tHV_from_R0(p, R0s)

out = DataFrame([Float64[],Float64[], Float64[], Float64[]], [:H, :V, :t, :R0])

for s in eachindex(THVs)
  R0 = R0s[s]
  p[2]=THVs[s] # Thv value
  # Define ODE problem (no stochasticity)
  prob = ODEProblem(f,u0, tspan, p)
  # Solve 
  sol = solve(prob)

  states = sol.u
  times = sol.t
  traj_det = DataFrame(sol)
  temp = reduce(vcat, transpose(states))
  temp2 = reduce(hcat, [temp, times])
  R0_val = repeat([R0], size(temp)[1])
  temp3 = reduce(hcat, [temp2, R0_val])
  temp4 = DataFrame(temp3, [:H, :V, :t, :R0])

  append!(out,temp4)

end

# Save trajectories
cd("$(homedir())/Documents/Github/dahlin-noisy-mbds/results") do
  CSV.write("deterministic_trajectories.csv", out, transform = (col,val) -> something(val, missing))
end


# Set up SDE problem_____________________________________

# Callback functions
# Callback 1: if infected host pop. goes negative, replace its value with zero
condition1(u, t, integrator) = u[1] < 0
function affect1!(integrator)
  integrator.u[1] = 0
  integrator.p[6] = 0
end
cb1 = ContinuousCallback(condition1, affect1!, save_positions=(false,true))
# Callback 2: if infected vector pop. goes negative, replace its value with zero
condition2(u, t, integrator) = u[2] < 0
function affect2!(integrator)
  integrator.u[2] = 0
  integrator.p[7] = 0
end
cb2 = ContinuousCallback(condition2, affect2!, save_positions=(false,true))
# Callback 4: if infected host pop. exceeds carrying capacity (10000), replace its value with carrying capacity
condition4(u, t, integrator) = u[1] > p[4]
function affect4!(integrator)
  integrator.u[1] = p[4]
  integrator.p[2] = 0
  # @show integrator.u[1], integrator.t 
end
cb4 = ContinuousCallback(condition4, affect4!, save_positions=(false,true))
# Callback 5:  if infected vector pop. exceeds carrying capacity (100000), replace its value with carrying capacity
condition5(u, t, integrator) = u[2] > p[5]
function affect5!(integrator)
  integrator.u[2] = p[5]
  integrator.p[3] = 0
end
cb5 = ContinuousCallback(condition5, affect5!, save_positions=(false,true))
# Combine callbacks
cbs = CallbackSet(cb1, cb2, cb4, cb5)

function condition(out, u, t, integrator)
  out[1] = u[1]
  out[2] = u[1] - p[4]
  out[3] = u[2]
  out[4] = u[2] - p[5]
end

function affect!(integrator, idx)
  if idx == 1
    # integrator.u[1] = 0
    integrator.p[6] = 0
  elseif idx == 2
    # integrator.u[1] = p[4]
    integrator.p[2] = 0
  elseif idx == 3
    # integrator.u[2] = 0
    integrator.p[7] = 0
  elseif idx == 4
    # integrator.u[2] = p[5]
    integrator.p[3] = 0
  end
end
# cbs = VectorContinuousCallback(condition, affect!, 4)


# Function: reduce to batches
function reduction(u,batch,I)
tmp = sum(cat(batch..., dims = 5), dims = 5)/length(I)
length(u) == 0 && return tmp, false
cat(u, tmp, dims = 5), false
end

# Define parameter space
# R0s = [1.05, 1.25, 2]
THVs = tHV_from_R0(p, R0s)
sigmas = [0.05, 0.25, 0.75, 0.95]
# Create grid to reduce number of for loops
const parm_grid = collect(Iterators.product(sigmas, R0s))

num_trajectories = 100

# Drift terms
function f(du,u,p,t)
  #(V,H) = u
  #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
  du[1] = (p[1] * p[2] * (p[4] - u[1]) / p[4]) * u[2] - p[6] * u[1]
  #du[1] = (b * τₕᵥ * (Nₕ - H) * V / Nₕ) - γₕ * H
  du[2] = (p[1] * p[3] * (p[5] - u[2]) / p[4]) * u[1] - p[7] * u[2]
  #du[2] = (b * τᵥₕ * (Nᵥ - V) * H / Nₕ) - μᵥ * V]
end
# Diffusion terms
function g(du,u,p,t)
  #(V,H) = u
  #(b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ, σ, αᵦ, αₜ, αₘ) = p
  du[1,1] = sqrt(max(0, (p[1] * p[2] * (p[4] - u[1])) / p[4]) * u[2] + p[6] * u[1])
  #du[1,1] = sqrt(max(0,(b * τₕᵥ * (Nₕ - H) * V / Nₕ) + γₕ * H))
  du[2,1] = 0
  du[1,2] = 0
  du[2,2] = sqrt(max(0, (p[1] * p[3] * (p[5] - u[2])) / p[4]) * u[1] + p[7] * u[2])
  #du[2,2] = sqrt(max(0,(b * τᵥₕ * (Nᵥ - V) * H / Nₕ) + μᵥ * V))
  du[1,3] = (p[1] * p[8] * p[9] * p[2] * (p[4] - u[1]) / p[4]) * u[2]
  #du[1,3] = σ * αᵦ * τₕᵥ * (Nₕ - H) * V / Nₕ
  du[2,3] = p[8] * ((p[1] * p[9] * p[3] * (p[5] - u[2]) / p[4]) * u[1] + (p[3] * p[10] * p[1] * (p[5] - u[2]) / p[4]) * u[1] + p[7] * p[11] * u[2])
  #du[2,3] = σ * (αᵦ * τᵥₕ * H * (Nᵥ - V) / Nₕ + αₜ * b * H * (Nᵥ - V) / Nₕ + αₘ * V)
end

# Diffusion parameters
noise_rate_prototype = [0.0 0.0 0.0; 0.0 0.0 0.0]

# Define SDE problem
out = DataFrame([Float64[],Float64[], Float64[], Float64[], Float64[], Float64[]],  [:H, :V, :t, :trajectory, :sigma, :R0])
for s in eachindex(parm_grid)
  
  p[8]=parm_grid[s][1] # sigma value
  R0=parm_grid[s][2] 
  tHV = tHV_from_R0(p, R0)
  p[2] = tHV# Thv value

  
  total_runs = num_trajectories
  #defines the SDE problem to be solved
  prob = SDEProblem(f, g, u0, tspan, p, noise_rate_prototype=noise_rate_prototype, callback = cbs)
  # Set up an ensemble problem-run 1000 iterations and save only the final time step values
  ensembleprob = EnsembleProblem(prob) 
  sol = DataFrame
  sol = DataFrame(solve(ensembleprob, EM(), callback = cbs, dt = dt, abstol = 1e-7, trajectories = total_runs, save_everystep = true, verbose = false, isoutofdomain = (u,p,t) -> any(x -> x < 0, u)))
  states = sol.u
  times = sol.t

  for i in 1:size(states)[1]
    temp = reduce(vcat, transpose(states[i]))
    temp2 = reduce(hcat, [temp, times[i]])
    trajectory_indexes = repeat([i], size(temp)[1])
    sigma_val = repeat([p[8]], size(temp)[1])
    R0_val = repeat([R0], size(temp)[1])
    temp3 = reduce(hcat, [temp2, trajectory_indexes, sigma_val, R0_val])
    temp4 = DataFrame(temp3, [:H, :V, :t, :trajectory, :sigma, :R0])

    append!(out,temp4)
  end
end


cd("$(homedir())/Documents/Github/dahlin-noisy-mbds/results") do
  CSV.write("illustrate_trajectories_proportionalSDE.csv", out, transform = (col,val) -> something(val, missing))
end