using Dates

model_type = "land" # "land" or "bucket"
const variable_list = ["swu"] # variable(s) you want to capture by adjusting your priors
const n_iterations = 10 # 1 iterations takes ~ 1.5 hour with current settings ((50, 15) resolution, 2 year simulation)
const spinup_period = Year(1)
const start_date = DateTime(2000, 12, 01) # this is the start of the forward model spinup
@assert month(start_date + spinup_period) == 12 "The start of your calibration period should be December."
const nelements = (101, 15) # resolution - (horizontal elements (lon, lat), vertical elements (soil discretization))
const dirname = "land_snow_zenith_hires" # ideally, describe your calibration in a few words
const caldir = joinpath("output", dirname) # you might want to save somewhere else than login
import ClimaLand
model_dir = joinpath(pkgdir(ClimaLand), "experiments", "calibration")

# Don't forget to adjust your priors and forward_model files.

# Most of the time users will just need to modify the settings above
# In future commits, more flexibility should be passed as argument above,
# for example, add Noise choice (inter-annual variance, flat per variable...), or simulation_period, ...

ENV["JULIA_WORKER_TIMEOUT"] = "1000.0" # This is a ClimaCalibrate setting. Wait time for workers.

using Dates
using Distributed
import EnsembleKalmanProcesses as EKP
import Random
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
rng_seed = 2
rng = Random.MersenneTwister(rng_seed)
FT = Float64

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments",
        "calibration",
        string("priors_", model_type, ".jl"),
    ),
)
# potentially we could add loss-parameters pairing https://clima.github.io/EnsembleKalmanProcesses.jl/dev/update_groups/

ekp_process = EKP.Unscented(prior)
ensemble_size = ekp_process.N_ens

# Config for workers
#addprocs(CAL.SlurmManager())
CAL.add_workers(ensemble_size; cluster = :auto, device = :gpu, time = 700)

@everywhere using Dates
@everywhere using ClimaLand
# Send variables from main process to all workers
@everywhere global model_type
@everywhere global nelements
@everywhere global caldir
@everywhere global start_date
@everywhere global model_dir
for pid in workers()
    @spawnat pid begin
        global model_type = Main.model_type
        global nelements = Main.nelements
        global caldir = Main.caldir
        global start_date = Main.start_date
        global model_dir = Main.model_dir
    end
end

@everywhere begin
    include(joinpath(model_dir, string("forward_model_", model_type, ".jl")))
end

# Locations used for calibration (currently all coordinates on land):
include(joinpath(model_dir, "make_training_locations.jl"))
training_locations = make_training_locations(nelements)
include(joinpath(model_dir, "observationseries_era5.jl"))

n_locations = length(training_locations)
n_variables = length(variable_list)
n_time_points = 4 # 4 seasons (and not, for example, 12 months)
l_obs = n_time_points * n_variables * n_locations

include(joinpath(model_dir, "observation_map.jl"))

# build observation from ClimaLand outputs - for all members
function CAL.observation_map(iteration)
    single_member_dims = (l_obs,)
    G_ensemble = Array{Float64}(undef, single_member_dims..., ensemble_size)

    for m in 1:ensemble_size
        member_path = path_to_ensemble_member(caldir, iteration, m)
        simdir_path =
            joinpath(member_path, "global_diagnostics", "output_active")
        simdir = SimDir(simdir_path)
        G_ensemble[:, m] .=
            process_member_data(simdir, training_locations, variable_list)
    end

    return G_ensemble
end

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments",
        "long_runs",
	"leaderboard",
	"leaderboard.jl"),
)
function CAL.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    plot_output_path = CAL.path_to_iteration(output_dir, iteration)

    member1_path = CAL.path_to_ensemble_member(output_dir, iteration, 1)
    diagnostics_folder_path =
        joinpath(member1_path, "global_diagnostics/output_0000")

    compute_monthly_leaderboard(
        plot_output_path,
        diagnostics_folder_path,
        "ERA5",
    )
    plot_constrained_params_and_errors(plot_output_path, ekp, prior)
end

function plot_constrained_params_and_errors(output_dir, ekp, prior)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
    for i in 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size + 1], ekp)
    EKP.Visualize.plot_error_over_time(fig[1, dim_size + 2], ekp)
    CairoMakie.save(
        joinpath(output_dir, "constrained_params_and_error.png"),
        fig,
    )
    return nothing
end
# Build the UTKI object - this is where you set EKP configurations
utki = EKP.EnsembleKalmanProcess(
    observationseries,
    EKP.TransformUnscented(prior, impose_prior = true);
    verbose = true,
    rng,
    scheduler = EKP.DataMisfitController(terminate_at = 100),
)

# Run the calibration via ClimaCalibrate using the arguments built above
CAL.calibrate(CAL.WorkerBackend, utki, n_iterations, prior, caldir)
