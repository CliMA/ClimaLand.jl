# The script included below contains the following functions:
# rand_locations - returns n random locations on land (lon, lat)
# setup_prob - config for the bucket model
# target_and_locations - returns locations on land and ERA5 LHF SHF obs at those locations
using ClimaLand
dir = pkgdir(ClimaLand)
# locations, observations = target_and_locations()
# Should observations be a vector or EKP.Observation type?

# Now, we need a forward_model for ClimaCalibrate.
# forward_model runs the model and generates a ClimaDiagnostic output directory
# that will be used to generate the observation_map.
# note that forward_model needs the setup_prob function defined above.
# IMPORTANT: forward_model needs to be defined as CAL.forward_model, as it adds a method
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
using Distributed
addprocs(CAL.SlurmManager())

@everywhere begin
    import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
    import ClimaCalibrate as CAL

    using ClimaLand
    caldir = "calibration_output"
    dir = pkgdir(ClimaLand)

    include(joinpath(dir,"experiments/calibration/global_bucket/climacalibrate_bucket/bucket_target_script.jl"))
    include(joinpath(dir,"experiments/calibration/global_bucket/climacalibrate_bucket/forward_model.jl"))
end

# Now we include the script observation_map.jl, which contains two functions:
# process_member_data - makes the parameter to data map we want from the forward_model output for one ensemble member
# observation_map - makes our G_ensemble, the parameter to data map for all n_ensemble parameters
# IMPORTANT: observation_map needs to be defined as CAL.observation_map, as it adds a method
include(joinpath(dir,"experiments/calibration/global_bucket/climacalibrate_bucket/observation_map.jl"))

# Now that our functions are defined, we still need some to specify some arguments they require.
# Parameters prior
# Note: this could be provided in a .toml file. Not sure if it has to.
prior_κ_soil = EKP.constrained_gaussian("κ_soil", 2, 1, 0, Inf);
prior_ρc_soil = EKP.constrained_gaussian("ρc_soil", 4e6, 2e6, 0, Inf);
prior_f_bucket = EKP.constrained_gaussian("f_bucket", 0.5, 0.3, 0, 1);
prior_W_f = EKP.constrained_gaussian("W_f", 0.4, 0.4, 0, Inf);
prior_p = EKP.constrained_gaussian("p", 2, 1, 1, Inf);
prior_z_0m = EKP.constrained_gaussian("z_0m", 0.01, 0.1, 0, Inf);
prior = EKP.combine_distributions([
    prior_κ_soil,
    prior_ρc_soil,
    prior_f_bucket,
    prior_W_f,
    prior_p,
    prior_z_0m,
]);

ensemble_size = 10
n_iterations = 5
noise = 1.0*EKP.I # Should work, but this should be covariance of each month from observation (ERA5)

# This folder will contains ClimaCalibrate outputs - parameters ensemble at each iterations, etc.
caldir = "calibration_output"

# Note: what is the best way to test this?
# Should this script be run as a slurm job? (I expect it takes ~ 5 hours or so to complete)
CAL.calibrate(
              CAL.WorkerBackend,
              ensemble_size,
              n_iterations,
              observations,
              noise,
              prior,
              caldir # to do: copy a tree of that folder structure
             )
