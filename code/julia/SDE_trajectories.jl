# Create example trajectories to illustrate how dynamics as noise type and magnitude are varied
include("RM_Reflected_SDE.jl")


# Sets of R0 - sigma
example_R0s = [0.95, 1.375, 2.05] # 1.5, 
example_sigmas = [0.05, 0.55, 1] # 0.3, 
example_Thvs = Thv_from_R0(q, example_R0s)

example_values = [(Thv, sigma) for Thv in example_Thvs, sigma in example_sigmas]

# Get solution trajectories representative of each set of values
num_sims = 1000::Int

trajectories_all_noise = run_sims(dF_det!, dF_stoch!, num_sims, example_values)
trajectories_all_noise[!,:type] = fill("All_noise", size(trajectories_all_noise)[1])

trajectories_no_demo = run_sims(dF_det_no_demo!, dF_stoch_no_demo!, num_sims, example_values)
trajectories_no_demo[!,:type] = fill("No_demographic", size(trajectories_no_demo)[1])



results = DataFrame(time = Float64[], Thv = Float64[], sigma = Float64[], run = Int[], H = Float64[], V = Float64[])
for i in ProgressBar(eachindex(example_values))

	## Set up ODE problem
	prob = ODEProblem(dF_det!, u0, timespan, example_values[i], callback = cbs)
	# ensembleprob = EnsembleProblem(prob)
	## Run SDE solver
	sol = solve(prob, Tsit5(); saveat = 30.0f0)

	# Collect results into a tidy DataFrame
	# for run_id in 1:num_runs
	trajectory = sol

	results_append = DataFrame(
		time = trajectory.t,
		Thv = fill(example_values[i][1], length(trajectory.t)),
		sigma = fill(example_values[i][2], length(trajectory.t)),
		run = fill(1, length(trajectory.t)),
		H = [u[1] for u in trajectory.u],
		V = [u[2] for u in trajectory.u],
	)

	# Append to the main results DataFrame
	append!(results, results_append)

end
trajectories_deterministic = results
trajectories_deterministic[!,:type] = fill("Deterministic", size(trajectories_deterministic)[1])

all_trajectories = append!(trajectories_deterministic, trajectories_all_noise, trajectories_no_demo)
CSV.write(joinpath(pwd(), "data", "comparison_trajectories.csv"), all_trajectories)