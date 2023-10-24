using DifferentialEquations
using Plots
using DataFrames
using Statistics
#using StatsBase
#using Suppressor # to suppress warnings from SDE solver
#using ProgressBars # to gauge how long the solver will take
#using TimerOutputs # to actually time which commands are taking the most time
using CSV
# using PlotlyJS


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
function tVH_from_R0(p, R0) 
  (b, τₕᵥ, τᵥₕ, Nₕ, Nᵥ, γₕ, μᵥ) = p
  tVH = (R0.^2) / ((b.^2 * τₕᵥ * Nᵥ) / (Nₕ * γₕ * μᵥ))
  return(tVH)
end

#Parameters: biting rate, THV, TVH, NH, NV, recovery, mortality, sigma, alphab, alpham, alphat
p = [0.3, 0.02450000, 0.5, 10000, 100000, 0.1, 0.1, 0.0125, 1.0, 1.0, 1.0]



#Rate of transmission from hosts to vectors-We varied this parameter to get different R0 values
#Code for calculating Tvh for a given R0 value is in the main R code
# TVHs=[0.01250000, 0.02005556, 0.02450000, 0.03472222, 0.08888889,  0.35555556, 0.93888889]#R0=0.75, 0.95, 1.05, 1.25, 2, 4, 6.5
R0s = [0.75, 0.95, 1.05, 1.25, 2, 4, 6.5]
TVHs = tVH_from_R0(p, R0s)
#[0.02005556, 0.02450000, 0.03472222, 0.08888889, 0.35555556, 0.93888889] # These are for R0=0.95, 1.05, 1.25, 2, 4, 6.5
#[ 0.02450000, 0.03472222, 0.05000000, 0.08888889, 0.35555556] #These are for R0=1.05, 1.25,1.5, 2, 4
#[0.01250000, 0.02005556, 0.02450000, 0.35555556] These are for R0=0.75, 0.95, 1.05, 4
#Now varying R0 by TVH since that is not a parameter we want to play with for sensitivity analysis

#initial conditions->10 infected mosquitos to star
u0 = [0.0, 10.0]
#10 year simulation with time step of 1 day
tspan = [0.0, 10*365]
dt = 1


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

  ##callback 4 makes Host population value 10000 if it exceeds this
function condition4(u, t, integrator)
      u[1]>10000
end
function affect4!(integrator)
      integrator.u[1] = 10000
end
cb4 = DiscreteCallback(condition4, affect4!)
  #callback 5 makes Vector population value 100000 if it exceeds this
function condition5(u, t, integrator)
    u[2]>100000
end
function affect5!(integrator)
    integrator.u[2] = 100000
end
cb5 = DiscreteCallback(condition5, affect5!)
cbs = CallbackSet(cb1, cb2, cb4, cb5)

  
prob = ODEProblem(f,u0, tspan, p, callback=cbs)
sol = solve(prob)
plot(sol, layout=(3,1), legend=true)
hline!([10])

p[3]=TVHs[2]
prob = ODEProblem(f,u0, tspan, p, callback=cbs)
sol = solve(prob)
plot(sol, layout=(3,1), legend=true)
hline!([10])

p[3]=TVHs[3]
prob = ODEProblem(f,u0, tspan, p, callback=cbs)
sol = solve(prob)
plot(sol, layout=(3,1), legend=true)
hline!([10])

p[3]=TVHs[4] #R0=1.25
prob = ODEProblem(f,u0, tspan, p, callback=cbs)
sol = solve(prob)
plot(sol, layout=(3,1), legend=true)
hline!([10])

traj_det = DataFrame(sol)
cd("/Users/karinebey/Documents/GitHub/dahlin-noisy-mbds/") do
  CSV.write("julia_trajectories_1.05deterministic.csv", traj_det, transform = (col,val) -> something(val, missing))
end

p[3]=TVHs[5]
prob = ODEProblem(f,u0, tspan, p, callback=cbs)
sol = solve(prob)
plot(sol, layout=(3,1), legend=true)
hline!([10])

p[3]=TVHs[6]
prob = ODEProblem(f,u0, tspan, p, callback=cbs)
sol = solve(prob)
plot(sol, layout=(3,1), legend=true)
hline!([10])

p[3]=TVHs[7]
prob = ODEProblem(f,u0, tspan, p, callback=cbs)
sol = solve(prob)
plot(sol, layout=(3,1), legend=true)
hline!([10])



#Stochastic trajectories R0=1.25; sigma=0.12 and 0.0

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

end

p = [0.3, 0.08888889, 0.5, 10000, 100000, 0.1, 0.1, 0.4, 1.0, 1.0, 1.0]
u0 = [0.0, 10.0]
tspan = [0.0, 10*365]
noise_rate_prototype = [0.0 0.0 0.0; 0.0 0.0 0.0]
dt = 1

#callbacks make value 0 if goes negative

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

##callback 4 makes Host population value 10000 if it exceeds this
function condition4(u, t, integrator)
    u[1]>10000
end
function affect4!(integrator)
    integrator.u[1] = 10000
end
cb4 = DiscreteCallback(condition4, affect4!)
#callback 5 makes Vector population value 100000 if it exceeds this
function condition5(u, t, integrator)
  u[2]>100000
end
function affect5!(integrator)
  integrator.u[2] = 100000
end
cb5 = DiscreteCallback(condition5, affect5!)
#This callback terminates the simulation if the outbreak dies out defined as less than 1 host and vector infected

cbs = CallbackSet(cb1, cb2, cb4, cb5)


prob = SDEProblem(f, g, u0, tspan, p, noise_rate_prototype=noise_rate_prototype, callback = cbs)#, callback = PositiveDomain())
#ensembleprob = EnsembleProblem(prob) 
#sol = solve(ensembleprob, EM(), dt = dt, EnsembleSplitThreads(), trajectories = 100, verbose = false)#, isoutofdomain = (u,p,t) -> any(x -> x < 0, u))
sol = solve(prob, EM(), dt = dt, verbose = false)

dftemp = empty

dftemp = DataFrame(sol)


cd("/Users/karinebey/Documents/GitHub/dahlin-noisy-mbds/") do
  CSV.write("traj_stoc_2_0.4_10.csv", dftemp, transform = (col,val) -> something(val, missing))
end


cd("/Users/karinebey/Documents/GitHub/NoisyMBDs") do
  CSV.write("traj_6.5envstoc_2.csv", dftemp, transform = (col,val) -> something(val, missing))
end
cd("/Users/karinebey/Documents/GitHub/NoisyMBDs") do
  CSV.write("traj_6.5envstoc_3.csv", dftemp, transform = (col,val) -> something(val, missing))
end
cd("/Users/karinebey/Documents/GitHub/NoisyMBDs") do
  CSV.write("traj_6.5envstoc_4.csv", dftemp, transform = (col,val) -> something(val, missing))
end
cd("/Users/karinebey/Documents/GitHub/NoisyMBDs") do
  CSV.write("traj_6.5envstoc_5.csv", dftemp, transform = (col,val) -> something(val, missing))
end
