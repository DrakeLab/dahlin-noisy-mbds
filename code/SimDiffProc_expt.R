library(Sim.DiffProc)

# Htemp <- Htemp + (b * tHV * Vtemp * (NH - Htemp) / NH - gH * Htemp) * Dt +
#   sqrt(b * tHV * Vtemp * (NH - Htemp) / NH + gH * Htemp) * W1inc +
#   sigma * (alphab * Vtemp * tHV * (NH - Htemp) / NH) * W3inc
# 
# Vtemp <- Vtemp + (b * tVH * Htemp * (NV - Vtemp) / NH - mV * Vtemp) * Dt +
#   sqrt(b * tVH * Htemp * (NV - Vtemp) / NH + mV * Vtemp) * W2inc +
#   sigma * (alphab * Htemp * tVH * (NV - Vtemp) / NH +
#              alphat * Htemp * b * (NV - Vtemp) / NH +
#              alpham * Vtemp) * W3inc

# Parameters -------------------------
b <- 0.3 / 4.5 # biting rate # KD: R0 equals either 1.05 or 4.74
tHV <- 0.5 # transmission prob. from V to H
tVH <- 0.5 # transmission prob. from H to V
NH <- 1000 # host population size
NV <- 10000 # vector population size
gH <- 0.1 # recovery rate of hosts
mV <- 0.1 # mortality of vectors
#** Environmental noise ----
sigma <- seq(0, 1, length.out = 21) # strength of env. noise
alphab <- 1 # biting env. noise on/off
alphat <- 1 # trans. prob. env. noise on/off
alpham <- 1 # mort. env. noise on/off

R0 <- sqrt((b^2 * tHV * tVH * NV) / (NH * gH * mV))

Hstar <- max((NH * (b^2 * NV * tHV * tVH - NH * mV * gH)) / (b * tVH * (b * NV * tHV + NH * gH)), 1000)
Vstar <- max((b^2 * NV * tVH * tHV - NH * mV * gH) / (b * tHV * (b * tVH + mV)), 10000)

# Initial conditions------------------------------------------------------------
x0 = 0; y0 = 100
tMax <- 4 * 365 # 500
dt <- 0.01
tspan <- seq(0.0, tMax, dt)
parms.vec <- c(b = b, tHV = tHV, tVH = tVH, NH = NH, NV = NV, gH = gH, mV = mV)

# Define SDEs-------------------------------------------------------------------
fx <- expression(b * tHV * y * (NH - x) / NH - gH * x, 
                 b * tVH * x * (NV - y) / NH - mV * y)

gx <- expression(sqrt(b * tHV * y * (NH - x) / NH + gH * x), 
                 sqrt(b * tVH * x * (NV - y) / NH + mV * y))

N = tMax/dt


set.seed(1234)

# Numerical solution of SDE-----------------------------------------------------
system.time(mod2d <- snssde2d(
  drift = fx, # drift term
  diffusion = gx, # diffusion term
  t0 = 0, # initial time
  T = tMax, # final time
  N = tMax/dt, # number of simulation steps (preserves dt step size)
  M = 20,  # number of trajectories to sample
  x0 = c(x0, y0), # initial values
  corr = NULL # correlation structure of Brownian motion
)
) 
# Plots-------------------------------------------------------------------------

# Trajectories over time
plot(mod2d)

# X vs Y
plot2d(mod2d, type = "n") 
points2d(mod2d, col = rgb(0, 100, 0, 50, maxColorValue = 255), pch = 16)

# Marginal density
denM <- dsde2d(mod2d, pdf="M", at = tMax) # Calculate marginal density at time endpoint
plot(denM, main="Marginal Density")

# Joint density
denJ <- dsde2d(mod2d, pdf = "J", n = 100, at = tMax/10)
plot(denJ, display = "contour", main = "Bivariate Transition Density at time t=tmax")

# Trying a different package...-------------------------------------------------

library(diffeqr)
library(JuliaCall)
julia <- julia_setup(JULIA_HOME = "C:\\Users\\kd99491\\AppData\\Local\\Programs\\Julia-1.9.1\\bin")
JuliaCall::julia_install_package_if_needed("DifferentialEquations")
JuliaCall::julia_library("DifferentialEquations")

# Set up differential equation
de <- diffeqr::diffeq_setup()

# fx <- expression(b * tHV * y * (NH - x) / NH - gH * x, 
#                  b * tVH * x * (NV - y) / NH - mV * y)
# 
# gx <- expression(sqrt(b * tHV * y * (NH - x) / NH + gH * x), 
#                  sqrt(b * tVH * x * (NV - y) / NH + mV * y))

# Parameter indexing:
# b = 1, tHV = 2, tVH = 3, NH = 4, NV = 5, gH = 6, mV = 7,
# sigma = 8, alpha_b = 9, alpha_tau = 10, alpha_mu = 11

f <- JuliaCall::julia_eval("
function f(du,u,p,t)
  du[1] = (p[1] * p[2] * (p[4] - u[1]) * u[2] / p[4]) - p[6] * u[1]
  du[2] = (p[1] * p[3] * (p[5] - u[2]) * u[1] / p[4]) - p[7] * u[2]
end")
g <- JuliaCall::julia_eval("
function g(du,u,p,t)
  du[1,1] = sqrt(max(0,(p[1] * p[2] * (p[4] - u[1]) * u[2] / p[4]) + p[6] * u[1]))
  du[2,1] = 0
  du[1,2] = 0
  du[2,2] = sqrt(max(0,(p[1] * p[3] * (p[5] - u[2]) * u[1] / p[4]) + p[7] * u[2]))
  du[1,3] = p[8] * p[9] * p[2] * (p[4] - u[1]) * u[2] / p[4]
  du[2,3] = p[8] * (p[9] * p[3] * (p[5] - u[2]) * u[1] / p[4] + p[10] * p[1] * (p[5] - u[2]) / p[4] + p[11] * u[2])
end")
p <- c(b, tHV, tVH, NH, NV, gH, mV,
       0.05, 1, 1, 1) 
u0 <- c(0.0, 100.0)
tspan <- c(0.0,10*365)
noise_rate_prototype <- matrix(c(0.0,0.0,0.0,0.0, 0.0, 0.0), nrow = 2, ncol = 3)
trajectories <- 10

JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("p", p)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("trajectories", trajectories)
JuliaCall::julia_assign("noise_rate_prototype", noise_rate_prototype)

prob <- JuliaCall::julia_eval("SDEProblem(f, g, u0, tspan, p, noise_rate_prototype=noise_rate_prototype)")
JuliaCall::julia_assign("prob", prob)

# ensembleprob <- JuliaCall::julia_eval("EnsembleProblem(prob, trajectories = 10)")
# JuliaCall::julia_assign("ensembleprob", ensembleprob)
# sol <- JuliaCall::julia_eval("solve(ensembleprob, EnsembleThreads(), trajectories = 10)")
# JuliaCall::julia_assign("sol", sol)
# summ <- JuliaCall::julia_eval("EnsembleSummary(sol, 0:0.01:1)")


sol <- de$solve(prob) # this is where the solver actually works. may take some time
udf <- as.data.frame(t(sapply(sol$u,identity)))

plotly::plot_ly(udf, x = ~V1, y = ~V2, type = 'scatter', mode = 'lines')




# f <- JuliaCall::julia_eval("
# function f(du,u,p,t)
#   du[1] = 10.0*(u[2]-u[1])
#   du[2] = u[1]*(28.0-u[3]) - u[2]
#   du[3] = u[1]*u[2] - (8/3)*u[3]
# end")
# g <- JuliaCall::julia_eval("
# function g(du,u,p,t)
#   du[1,1] = 0.3u[1]
#   du[2,1] = 0.6u[1]
#   du[3,1] = 0.2u[1]
#   du[1,2] = 1.2u[2]
#   du[2,2] = 0.2u[2]
#   du[3,2] = 0.3u[2]
# end")
# u0 <- c(1.0,0.0,0.0)
# tspan <- c(0.0,100.0)
# noise_rate_prototype <- matrix(c(0.0,0.0,0.0,0.0,0.0,0.0), nrow = 3, ncol = 2)
# 
# JuliaCall::julia_assign("u0", u0)
# JuliaCall::julia_assign("p", p)
# JuliaCall::julia_assign("tspan", tspan)
# JuliaCall::julia_assign("noise_rate_prototype", noise_rate_prototype)
# prob <- JuliaCall::julia_eval("SDEProblem(f, g, u0, tspan, p, noise_rate_prototype=noise_rate_prototype)")
# 
# sol <- de$solve(prob) # takes a long time
# udf <- as.data.frame(t(sapply(sol$u,identity)))
# plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')






