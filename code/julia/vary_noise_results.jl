#### Get results and save ####
include("RM_Reflected_SDE.jl")
using CodecZlib

test = collect_outputs(dF_det!, dF_stoch!, 10, parameter_values[1:2])
test2 = raw_outputs(dF_det!, dF_stoch!, 10, parameter_values[1:2])
# using Profile # ProfileView
# @profile collect_all_no_demo = collect_outputs(dF_det_no_demo!, dF_stoch_no_demo!, 20, parameter_values)
# test_no_demo = collect_outputs(dF_det_no_demo!, dF_stoch_no_demo!, 1, parameter_values[1:2])
test_det = collect_outputs_det(dF_det!, Thvs[1])

num_sims = 100_000::Int #10_000::Int# number of simulations to Run

function save_gzip(data, filename)
    open(joinpath(dirname(dirname(pwd())), "data", filename), "w") do io
        gzip_io = GzipCompressorStream(io)
        CSV.write(gzip_io, data, append = true, writeheader = true)
        close(gzip_io)
    end
end

using CSV
using CodecZlib
using DataFrames

function append_and_save_gzip(data::DataFrame, filename::String)
    filepath = joinpath(dirname(dirname(pwd())), "data", filename)

    # Load existing data if the file exists
    if isfile(filepath)
        open(filepath) do io
            gzip_io = GzipDecompressorStream(io)
            existing_data = CSV.read(gzip_io, DataFrame)
            close(gzip_io)
            data = vcat(existing_data, data)
        end
    end

    # Save the combined data
    open(filepath, "w") do io
        gzip_io = GzipCompressorStream(io)
        CSV.write(gzip_io, data)
        close(gzip_io)
    end
end



Thvs_for_trajectories = Thv_from_R0(q, [0.95, 1.05, 1.25, 1.375, 2, 3, 4, 4.625]) # used to vary R0
sigmas_for_trajectories = [0, 0.1, 0.25, 0.5, 1, 1.5]
parameter_values_for_trajectories = [(Thv, sigma) for Thv in Thvs_for_trajectories, sigma in sigmas_for_trajectories]


## Simulations without any noise ##
print("Simulations with no noise at all")
collect_all_det = collect_outputs_det(dF_det!, R0s)
save_gzip(collect_all_det, "collect_all_outputs_det.csv.gz")
finalize(collect_all_det)
collect_all_det = nothing
GC.gc()

parameter_values_det = [(Thv, 0) for Thv in Thvs_for_trajectories]
trajectories_for_grid_plot_det = run_sims(dF_det_no_demo!, dF_stoch_no_demo!, 1::Int, parameter_values_det)
save_gzip(trajectories_for_grid_plot_det, "trajectories_for_grid_plot_det.csv.gz")
finalize(trajectories_for_grid_plot_det)
GC.gc()


## Simulations without demographic noise ##
print("Simulations with no demographic noise")
collect_all_no_demo = collect_outputs(dF_det_no_demo!, dF_stoch_no_demo!, num_sims, parameter_values)
save_gzip(collect_all_no_demo, "collect_all_outputs_no_demo.csv.gz")
finalize(collect_all_no_demo)
collect_all_no_demo = nothing
GC.gc()

trajectories_for_grid_plot_no_demo = run_sims(dF_det_no_demo!, dF_stoch_no_demo!, 100::Int, parameter_values_for_trajectories)
save_gzip(trajectories_for_grid_plot_no_demo, "trajectories_for_grid_plot_no_demo.csv.gz")
finalize(trajectories_for_grid_plot_no_demo)
GC.gc()

# Get data for duration vs intensity scatter plots
Thvs_dur_peak = Thv_from_R0(q, 0.125f0:0.125f0:5f0) # used to vary R0
sigmas_dur_peak = 0.25f0:0.25f0:1.5f0
dur_peak_par_vals = [(Thv, sigma) for Thv in Thvs_dur_peak, sigma in sigmas_dur_peak]
duration_v_peak = raw_outputs(dF_det_no_demo!, dF_stoch_no_demo!, 1_000, dur_peak_par_vals)
save_gzip(duration_v_peak, "dur_peak_no_demo.csv.gz")
collect_all = duration_v_peak
GC.gc()


## Simulations with all types of noise ##
println("Simulations with all noise")
chunk_size = 41
# Break up parameter values into smaller chunks
for i in 1:chunk_size:length(parameter_values)

    println("Simulating chunk ", i, " to ", i+chunk_size-1)

    chunk = parameter_values[i:minimum([i+chunk_size-1, length(parameter_values)])]
    # Get data for each set of parameter values
    collect_all = collect_outputs(dF_det!, dF_stoch!, num_sims, chunk)
    # Load previous file

    # Add rows to previous file

    # Save file
    append_and_save_gzip(collect_all, "collect_all_outputs.csv.gz")
    # save_gzip(collect_all, "collect_all_outputs.csv.gz")
    finalize(collect_all)
    collect_all = nothing
    GC.gc()
end

# CSV.write(joinpath(dirname(dirname(pwd())), "data", "collect_all_outputs.csv"), collect_all)


trajectories_for_grid_plot = run_sims(dF_det!, dF_stoch!, 1_000::Int, parameter_values_for_trajectories)
save_gzip(trajectories_for_grid_plot, "trajectories_for_grid_plot.csv.gz")
finalize(trajectories_for_grid_plot)
trajectories_for_grid_plot = nothing
GC.gc()
# CSV.write(joinpath(dirname(dirname(pwd())), "data", "trajectories_for_grid_plot.csv"), trajectories_for_grid_plot)

# Get data for duration vs intensity scatter plots
duration_v_peak = raw_outputs(dF_det!, dF_stoch!, 10_000, dur_peak_par_vals)
save_gzip(duration_v_peak, "dur_peak.csv.gz")
trajectories_for_grid_plot = duration_v_peak
GC.gc()