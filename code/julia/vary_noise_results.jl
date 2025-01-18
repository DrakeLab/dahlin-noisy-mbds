#### Get results and save ####
include("RM_Reflected_SDE.jl")

test = collect_outputs(dF_det!, dF_stoch!, 1, parameter_values)
test_no_demo = collect_outputs(dF_det_no_demo!, dF_stoch_no_demo!, 1, parameter_values)
test_det = collect_outputs_det(dF_det!, Thvs)

num_sims = 10000::Int# number of simulations to Run

# Simulations with all types of noise
print("Simulations with all noise")
collect_all = collect_outputs(dF_det!, dF_stoch!, num_sims, parameter_values)
CSV.write(joinpath(dirname(dirname(pwd())), "data", "collect_all_outputs.csv"), collect_all)

trajectories_for_grid_plot = run_sims(dF_det!, dF_stoch!, 50::Int, parameter_values)
CSV.write(joinpath(dirname(dirname(pwd())), "data", "trajectories_for_grid_plot.csv"), trajectories_for_grid_plot)

# Simulations without demographic noise
print("Simulations with no demographic noise")
collect_all_no_demo = collect_outputs(dF_det_no_demo!, dF_stoch_no_demo!, num_sims, parameter_values)
CSV.write(joinpath(dirname(dirname(pwd())), "data", "collect_all_outputs_no_demo.csv"), collect_all_no_demo)

trajectories_for_grid_plot_no_demo = run_sims(dF_det_no_demo!, dF_stoch_no_demo!, 50::Int, parameter_values)
CSV.write(joinpath(dirname(dirname(pwd())), "data", "trajectories_for_grid_plot_no_demo.csv"), trajectories_for_grid_plot_no_demo)

# Simulations without any noise
print("Simulations with no noise at all")
collect_all_det = collect_outputs_det(dF_det!, R0s)
CSV.write(joinpath(dirname(dirname(pwd())), "data", "collect_all_outputs_det.csv"), collect_all_det)

parameter_values_det = [(Thv, 0) for Thv in Thvs]
trajectories_for_grid_plot_det = run_sims(dF_det_no_demo!, dF_stoch_no_demo!, 1::Int, parameter_values_det)
CSV.write(joinpath(dirname(dirname(pwd())), "data", "trajectories_for_grid_plot_det.csv"), trajectories_for_grid_plot_det)
