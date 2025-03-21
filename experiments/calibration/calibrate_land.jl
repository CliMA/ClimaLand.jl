######## CONFIGS ########
# Most of the time users will just need to modify the settings below
# Our goal is to make it flexible and easy to run different calibration setups
# And to potentially iterate fast (e.g., low resolution, specific regions, few params, short simulation time...)
#########################

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

# Calibrate the land model with:
const variable_list = ["swu"] # variables you want to capture by adjusting your priors
const n_iterations = 10 # 1 iterations takes ~ 1.5 hour with current settings (resolution, 2 year simulation)
const spinup_period = Year(1)
# potentially we could add time_for_calibration (currently 1 year)

# Using the following priors:
include(joinpath(@__DIR__, "priors.jl"))
# potentially we could add loss-parameters pairing https://clima.github.io/EnsembleKalmanProcesses.jl/dev/update_groups/

# With the forward model on GPUS. Note: the forward model needs to be adjusted to read priors!
ekp_process = EKP.Unscented(prior)
ensemble_size = ekp_process.N_ens
# Config for workers
addprocs(
    CAL.PBSManager(ensemble_size), # simulations in parallel. User may change this to use less GPUs.
    q = "main",
    A = "UCIT0011",
    l_select = "1:ngpus=1:ncpus=4",
    l_walltime = "11:30:00",
    l_job_priority = "premium",
)
# ClimaLand (forward model) - needs to read priors!
@everywhere begin
    using Dates # needs to be called again for workers
    const nelements = (50, 10) # resolution - (horizontal elements (lon,lat), vertical elements (soil discretization))
    const start_date = DateTime(2008, 12, 01) # this is the start of the forward model spinup
    const caldir = "calibration_output_utki"
    using ClimaLand
    dir = pkgdir(ClimaLand)
    include(joinpath(dir, "experiments/calibration/forward_model.jl"))
end
@assert month(start_date + spinup_period) == 12 "The start of your calibration period should be December."

# And using those locations (currently all coordinates on land):
include(joinpath(@__DIR__, "make_training_locations.jl"))
training_locations = make_training_locations(nelements)
# potentially we can add regional runs or specific lon lat bands or filter (e.g., regions with snow)

# NOTE1: Don't forget to modify your forward model (land or bucket) to read and use your priors correctly.
# NOTE2: The noise is set in observationseries_era5.jl - adjust if needed.
# ^ current noise options: era5 inter-annual variance, era5 seasonal mean * factor, flat noise, weigh by lats.
# NOTE3: Currently everything is set to use seasonal averages. We could add option to use e.g., monthly or other.
# NOTE4: We could add option to calibrate single sites (change forward model, observationseries, observation_map).
# ^ maybe Julia and Thanhthanh SURF?



##########################
# Most of the time you won't need to change the code below.
# Unless you change EKP configurations, or other specific settings.
##########################


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
