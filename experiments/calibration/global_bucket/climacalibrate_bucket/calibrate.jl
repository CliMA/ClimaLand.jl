# The script included below contains the following functions:
# rand_locations - returns n random locations on land (lon, lat)
# setup_prob - config for the bucket model
# it also returns the target used for calibration, full_obs_era5, which is
# latent and sensible heat at the rand_locations stacked as a vector
include("bucket_target_script.jl")
observations = full_obs_era5

# Now, we need a forward_model for ClimaCalibrate.
# forward_model runs the model and generates a ClimaDiagnostic output directory
# that will be used to generate the observation_map.
# note that forward_model needs the setup_prob function defined above.
import ClimaCalibrate: forward_model, parameter_path
import ClimaCalibrate as CAL
include("forward_model.jl")

# Now we include the script observation_map.jl, which contains two functions:
# process_member_data - makes the loss function we want from the forward_model output
# observation_map - makes our G_ensemble, the loss for all n_ensemble parameters
include("observation_map.jl")

# Now that our functions are defined, we still need some to specify some arguments they require.
# Parameters prior
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
noise = EKP.I # not sure what to do here, I had noise embeded in my observation_map function

output_dir = "calibration_output" # Should I create this folder?

CAL.calibrate(
              CAL.WorkerBackend, # not sure what to do here
              ensemble_size,
              n_iterations,
              observations,
              noise,
              prior,
              caldir # to do: copy a tree of that folder structure
             )
