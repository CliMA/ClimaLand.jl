using Dates

model_type = "land" # "land" or "bucket"
const variable_list = ["swu"] # variable(s) you want to capture by adjusting your priors
const n_iterations = 10 # 1 iterations takes ~ 1.5 hour with current settings ((50, 15) resolution, 2 year simulation)
const spinup_period = Year(1)
const start_date = DateTime(2008, 12, 01) # this is the start of the forward model spinup
@assert month(start_date + spinup_period) == 12 "The start of your calibration period should be December."
const nelements = (50, 15) # resolution - (horizontal elements (lon, lat), vertical elements (soil discretization))
const dirname = "land_swu_zenith_only" # ideally, describe your calibration in a few words
const caldir = joinpath("/glade/campaign/univ/ucit0011/alexis/", dirname) # you might want to save somewhere else than login
# Don't forget to adjust your priors and forward_model files.



# Most of the time users will just need to modify the settings above
# In future commits, more flexibility should be passed as argument above,
# for example, add Noise choice (inter-annual variance, flat per variable...), or simulation_period, ...

ENV["JULIA_WORKER_TIMEOUT"] = "1000.0" # This is a ClimaCalibrate setting. Wait time for workers.

using Dates
using Distributed
import EnsembleKalmanProcesses as EKP
import Random
using ClimaLand
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
rng_seed = 2
rng = Random.MersenneTwister(rng_seed)
FT = Float64

include(joinpath(@__DIR__, string("priors_", model_type, ".jl")))
# potentially we could add loss-parameters pairing https://clima.github.io/EnsembleKalmanProcesses.jl/dev/update_groups/

ekp_process = EKP.Unscented(prior)
ensemble_size = ekp_process.N_ens

# Config for workers
CAL.add_workers(
    ensemble_size;
    cluster = :auto,
    device = :gpu,
    l_walltime = "11:30:00", # this only works for PBS. ClimaCalibrate will soon have a new more generalizable keyword for it to work on Slurm and PBS.
    l_job_priority = "premium", # this only works on PBS
)

@everywhere using Dates
@everywhere using ClimaLand
# Send variables from main process to all workers
@everywhere global model_type
@everywhere global nelements
@everywhere global caldir
@everywhere global start_date
for pid in workers()
    @spawnat pid begin
        global model_type = Main.model_type
        global nelements = Main.nelements
        global caldir = Main.caldir
        global start_date = Main.start_date
    end
end
@everywhere begin
    dir = pkgdir(ClimaLand)
    include(
        joinpath(
            dir,
            string("experiments/calibration/forward_model_", model_type, ".jl"),
        ),
    )
end

# Locations used for calibration (currently all coordinates on land):
include(joinpath(@__DIR__, "make_training_locations.jl"))
training_locations = make_training_locations(nelements)
# potentially we can add regional runs or specific lon lat bands or filter (e.g., regions with snow)

# NOTE: The noise is set in observationseries_era5.jl - adjust if needed.
# ^ current noise options: era5 inter-annual variance, era5 seasonal mean * factor, flat noise, weigh by lats.
# NOTE: Currently everything is set to use seasonal averages. We could add option to use e.g., monthly or other.
# NOTE: We could add option to calibrate single sites (change forward model, observationseries, observation_map).
# ^ maybe Julia and Thanhthanh SURF?

# observationseries - era5 data and noise object to compare to model output in EKP (to minimize the loss)
include(joinpath(@__DIR__, "observationseries_era5.jl"))

# l_obs is the length of Observation objects (observationseries for era5, observationmap for ClimaLand)
n_locations = length(training_locations)
n_variables = length(variable_list)
n_time_points = 4 # 4 seasons (and not, for example, 12 months)
l_obs = n_time_points * n_variables * n_locations

# build observation from ClimaLand outputs - for one member
include(joinpath(@__DIR__, "observation_map.jl"))

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
