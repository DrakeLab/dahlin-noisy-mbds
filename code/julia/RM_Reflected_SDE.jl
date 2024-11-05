# Numerical simulation of a Ross-Macdonald reflected SDE

# Load libraries
using DifferentialEquations
using Plots
using DataFrames
using DifferentialEquations.EnsembleAnalysis

# [] Later set this up as a function to parameters can be varied

# Set simulation parameters 
const num_sims = 1e5 # number of simulations to Run
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

const R0s = [6, 5, 4, 3, 2, 1.25, 1.2, 1.15, 1.1, 1.05, 1, 0.95, 0.75, 0.5, 0] # values of R0 to consider

function Thv_from_R0(q, R0) 
    (b, Tvh, Nh, Nv, muv, gammah) = q
    Thv = (R0.^2) / ((b.^2 * Tvh * Nv) / (Nh * gammah * muv))
    return(Thv)
  end

Thvs = Thv_from_R0(q,R0s) # used to vary R0

const sigmas = 0:0.05:1 # levels of environmental noise

# Initialize dataframes, if necessary

# Define the Ross-Macdonald Reflected SDE 

p = [b, Thvs[1], Tvh, Nh, Nv, muv, gammah,  sigmas[10], alphab, alphat, alpham] # !!! temp for testing. will vary from a dataframe later

## Deterministic part
function dF_det(du,u,p,t)
    # Name the variables
    (H,V) = u
    (b, Thv, Tvh, Nh, Nv, muv, gammah, sigma, alphab, alphat, alpham) = p
    
    du[1] = @views (b * Thv * (Nh - H) / Nh) * V - gammah * H
    du[2] = @views b * (Tvh / Nh) * H * (Nv - V) - muv * V
    return(du)
  end

## Stochastic part
function dF_stoch(du,u,p,t)
    # Name the variables
    (H,V) = u
    (b, Thv, Tvh, Nh, Nv, muv, gammah, sigma, alphab, alphat, alpham) = p
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
    if integrator.u[1] < 0.0
        integrator.u[1] = 0.0
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
    if integrator.u[2] < 0.0
        integrator.u[2] = 0.0
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
    integrator.u[2] = 0
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
    u[2] < 1 && u[1] < 1
end
affect_terminate!(integrator) = terminate!(integrator)
cb_terminate = ContinuousCallback(condition_terminate, affect_terminate!)

## Collect callbacks
# cbs = CallbackSet(cb_H0, cb_HNh, cb_V0, cb_VNv, cb_terminate)
cbs = CallbackSet(cb_H, cb_V, cb_terminate)

# Run simulations

## Set up ensemble SDE problem
prob = SDEProblem(dF_det, dF_stoch, u0, timespan, p, noise_rate_prototype = noise_rate_prototype, callback = cbs)
ensembleprob = EnsembleProblem(prob)


## Run SDE solver
sol = solve(prob, EM(), dt = 1) # works
# sol = solve(prob, SRA()) # steps too big, makes sqrt stuff go negative
# sol = solve(prob, SOSRA(), maxiters = 1e10) # needs high maxiters and still flat
# sol = solve(prob, SOSRA2()) # flat
# sol = solve(prob, SROCKEM(), dt = 1) # jumps over callbacks, so sqrt goes negative
# sol = solve(prob, ISSEM()) # jumps over callbacks, so sqrt goes negative
# sol = solve(prob, ISSEulerHeun())
sols = solve(ensembleprob, EM(), dt = 1, trajectories = 1)
sols = solve(ensembleprob, EM(), dt = deltat, EnsembleSplitThreads(), trajectories = 10000)
## Calculate summary statistics over time? [avoid this because it'll likely be much easier to do in R]


# Save simulations