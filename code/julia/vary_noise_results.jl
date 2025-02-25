#### Get results and save ####
include("RM_Reflected_SDE.jl")
using CodecZlib

test = collect_outputs(dF_det!, dF_stoch!, 1, parameter_values)
test_no_demo = collect_outputs(dF_det_no_demo!, dF_stoch_no_demo!, 1, parameter_values)
test_det = collect_outputs_det(dF_det!, Thvs)

num_sims = 10_000::Int# number of simulations to Run

function save_gzip(data, filename)
    open(joinpath(dirname(dirname(pwd())), "data", filename), "w") do io
        gzip_io = GzipCompressorStream(io)
        CSV.write(gzip_io, data)
        close(gzip_io)
    end
end

## Simulations without any noise ##
print("Simulations with no noise at all")
collect_all_det = collect_outputs_det(dF_det!, R0s)
save_gzip(collect_all_det, "collect_all_outputs_det.csv.gz")
finalize(collect_all_det)
collect_all_det = nothing
GC.gc()
# CSV.write(joinpath(dirname(dirname(pwd())), "data", "collect_all_outputs_det.csv"), collect_all_det)

parameter_values_det = [(Thv, 0) for Thv in Thvs]
trajectories_for_grid_plot_det = run_sims(dF_det_no_demo!, dF_stoch_no_demo!, 1::Int, parameter_values_det)
save_gzip(trajectories_for_grid_plot_det, "trajectories_for_grid_plot_det.csv.gz")
finalize(trajectories_for_grid_plot_det)
trajectories_for_grid_plot_det = nothing
GC.gc()
# CSV.write(joinpath(dirname(dirname(pwd())), "data", "trajectories_for_grid_plot_det.csv"), trajectories_for_grid_plot_det)


## Simulations without demographic noise ##
print("Simulations with no demographic noise")
collect_all_no_demo = collect_outputs(dF_det_no_demo!, dF_stoch_no_demo!, num_sims, parameter_values)
save_gzip(collect_all_no_demo, "t.csv.gz")
finalize(collect_all_no_demo)
collect_all_no_demo = nothing
GC.gc()
# CSV.write(joinpath(dirname(dirname(pwd())), "data", "collect_all_outputs_no_demo.csv"), collect_all_no_demo)

trajectories_for_grid_plot_no_demo = run_sims(dF_det_no_demo!, dF_stoch_no_demo!, 50::Int, parameter_values)
save_gzip(trajectories_for_grid_plot_no_demo, "trajectories_for_grid_plot_no_demo.csv.gz")
finalize(trajectories_for_grid_plot_no_demo)
trajectories_for_grid_plot_no_demo = nothing
GC.gc()
# CSV.write(joinpath(dirname(dirname(pwd())), "data", "trajectories_for_grid_plot_no_demo.csv"), trajectories_for_grid_plot_no_demo)


## Simulations with all types of noise ##
print("Simulations with all noise")
collect_all = collect_outputs(dF_det!, dF_stoch!, num_sims, parameter_values)
save_gzip(collect_all, "collect_all_outputs.csv.gz")
finalize(collect_all)
collect_all = nothing
GC.gc()
# CSV.write(joinpath(dirname(dirname(pwd())), "data", "collect_all_outputs.csv"), collect_all)

trajectories_for_grid_plot = run_sims(dF_det!, dF_stoch!, 50::Int, parameter_values)
save_gzip(trajectories_for_grid_plot, "trajectories_for_grid_plot.csv.gz")
finalize(trajectories_for_grid_plot)
trajectories_for_grid_plot = nothing
GC.gc()
# CSV.write(joinpath(dirname(dirname(pwd())), "data", "trajectories_for_grid_plot.csv"), trajectories_for_grid_plot)
