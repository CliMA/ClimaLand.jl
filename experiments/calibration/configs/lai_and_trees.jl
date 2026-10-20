# Calibration configuration for the optimal-LAI model against MODIS LAI and the
# CLM tree share
#
# This file is included by run_calibration.jl and defines:
#   - CALIBRATE_CONFIG (CalibrateConfig)
#   - NOISE_SCALARS (per-variable covariance scaling)
#   - get_calibration_prior() (combined prior distribution)
#
# To use this config, set the CALIBRATION_CONFIG env var when invoking
# run_calibration.sh, e.g.
#   CALIBRATION_CONFIG=lai_and_trees.jl bash experiments/calibration/run_calibration.sh
#
# The land model computes LAI with the optimal-LAI model (`ZhouOptimalLAIModel`) and
# its climate tree share, as the optimal-LAI long run does, from the spun-up
# optimal-LAI state of March 2008. Each member runs three years of high-resolution
# ERA5 (March 2008 to March 2011), and the third year is compared, as seasonal
# means, with MODIS LAI (`lai`) and with the tree share of the natural vegetation in
# the CLM5 surface data (`ftr`), both over natural vegetation only (cells with at
# most half of their land in crops), as in the vegetation leaderboard.
#
# The parameters were first calibrated in the offline emulator of the model
# (experiments/calibration/optimal_lai_emulator), whose values are the prior means.

"""Noise scalars for the covariance matrix of each observed variable.

These multiply the identity in `ScalarCovariance`, so they are variances.
- `lai`: (m^2 m^-2)^2. 0.5 is a standard deviation of about 0.7, the scale of the
  seasonal LAI errors.
- `ftr`: fraction^2. 0.2 is a standard deviation of about 0.45, twice the current
  error of the climate tree share (about 0.22), so that the tree share keeps its
  meaning without competing with LAI: its misfit is about a fifth of that of LAI.
"""
const NOISE_SCALARS = Dict("lai" => 0.5, "ftr" => 0.2)

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters.

Calibrates 9 parameters of the optimal-LAI model. Ensemble size for
TransformUnscented: 9 * 2 + 1 = 19 members.

- `optimal_lai_z_tree`, `optimal_lai_z_grass`: the leaf costs of trees and grasses
  (mol m^-2 yr^-1), which set the energy limit of LAI_max.
- `optimal_lai_maintenance_share`: the share of the leaf cost that scales with the
  days and temperature above freezing.
- `optimal_lai_sigma`: the departure from square-wave LAI dynamics.
- `optimal_lai_tree_retention`: the fraction of LAI_max that trees keep through
  the unfavourable season.
- `optimal_lai_tree_b0`, `_b_lai`, `_b_dry`, `_b_temp`: the logistic of the
  climate tree share. The priors are wide enough to reach the fit of this logistic
  to the CLM tree share on the model's 1° climate (-3.1, 0.34, -0.28, 0.12).

Left fixed: `optimal_lai_k`, degenerate with the leaf costs in the energy limit
1 - z/(k A0); `optimal_lai_alpha`, which earlier LAI calibrations settled at its
value; and the f0 curve and `optimal_lai_maintenance_q10`, which the emulator
could not identify.
"""
function get_calibration_prior()
    priors = [
        EKP.constrained_gaussian("optimal_lai_z_tree", 9.92, 3.0, 2.0, 40.0),
        EKP.constrained_gaussian(
            "optimal_lai_z_grass",
            154.0,
            50.0,
            30.0,
            500.0,
        ),
        EKP.constrained_gaussian(
            "optimal_lai_maintenance_share",
            0.43,
            0.15,
            0.0,
            1.0,
        ),
        EKP.constrained_gaussian("optimal_lai_sigma", 1.09, 0.25, 0.3, 3.0),
        EKP.constrained_gaussian(
            "optimal_lai_tree_retention",
            0.5,
            0.15,
            0.0,
            1.0,
        ),
        EKP.constrained_gaussian("optimal_lai_tree_b0", -1.07, 1.5, -8.0, 4.0),
        EKP.constrained_gaussian(
            "optimal_lai_tree_b_lai",
            0.0254,
            0.2,
            -0.5,
            1.5,
        ),
        EKP.constrained_gaussian(
            "optimal_lai_tree_b_dry",
            -0.372,
            0.15,
            -1.5,
            0.2,
        ),
        EKP.constrained_gaussian(
            "optimal_lai_tree_b_temp",
            0.0917,
            0.05,
            -0.1,
            0.4,
        ),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lai", "ftr"],
    minibatch_size = 1,
    n_iterations = 10,
    # One sample, the third year: MAM, JJA, SON 2010 and DJF 2010-2011 (Dec 1 is
    # the start of DJF), with extend = Month(3) running to March 1 2011.
    sample_date_ranges = [("2010-3-1", "2010-12-1")],
    extend = Dates.Month(3),
    # From March 1 2008, the date of the spun-up optimal-LAI state; the climate
    # totals that start empty (the moist season, the snow store) settle in a year.
    spinup = Dates.Year(2),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_lai_and_trees.jld2",
    model_type = ClimaLand.LandModel,
    prognostic_lai = true,
)
