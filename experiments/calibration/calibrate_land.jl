ENV["JULIA_WORKER_TIMEOUT"] = "1000.0"

using ClimaLand
dir = pkgdir(ClimaLand)
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
FT = Float64
using Distributed
using Random
rng_seed = 2
rng = Random.MersenneTwister(rng_seed)

addprocs(
    CAL.PBSManager(19), # n simulations in parallel
    q = "main", # not prio
    A = "UCIT0011",
    l_select = "1:ngpus=1:ncpus=4",
    l_walltime = "11:30:00",
    l_job_priority = "premium",
)

@everywhere begin
    import ClimaCalibrate:
        forward_model, parameter_path, path_to_ensemble_member
    import ClimaCalibrate as CAL

    using ClimaLand
    caldir = "calibration_output_utki"
    dir = pkgdir(ClimaLand)
    nelements = (50, 10) # resolution for model and era5

    include(joinpath(dir, "experiments/calibration/forward_model.jl"))
end

import EnsembleKalmanProcesses as EKP
include(joinpath(dir, "experiments", "calibration", "priors.jl"))

# ensemble_size = 40 # ideally 40
n_iterations = 10 # try 5, ideally 10
nelements = (50, 10) # resolution for model and era5

# training_locations
include(
    joinpath(dir, "experiments", "calibration", "make_training_locations.jl"),
)
path_to_grid = "calibration_output_utki_sample/iteration_000/member_001/global_diagnostics/output_active/"
training_locations = make_training_locations(path_to_grid)

# observationseries
include(joinpath(dir, "experiments/calibration/observationseries_era5.jl"))

n_locations = length(training_locations)
n_variables = 4 # LHF, SHF, SWU, LWU
n_time_points = 4 # 4 seasons (and not, for example, 12 months)
l_obs = n_time_points * n_variables * n_locations
ekp_process = EKP.Unscented(prior)
ensemble_size = ekp_process.N_ens

caldir = "calibration_output_utki"

include(joinpath(dir, "experiments/calibration/observation_map.jl"))

utki = EKP.EnsembleKalmanProcess(
    # EKP.construct_initial_ensemble(rng, prior, ensemble_size),
    observationseries,
    # EKP.Inversion();
    EKP.TransformUnscented(prior, impose_prior = true);
    verbose = true,
    # localization_method=EKP.Localizers.NoLocalization(),
    rng = rng,
    # scheduler=EKP.DefaultScheduler(1.0)
    scheduler = EKP.DataMisfitController(terminate_at = 100),
)

CAL.calibrate(CAL.WorkerBackend, utki, n_iterations, prior, caldir)
