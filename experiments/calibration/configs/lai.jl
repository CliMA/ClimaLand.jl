# Calibration config: optimal-LAI parameters against MODIS LAI.
# Included by run_calibration.jl; defines CALIBRATE_CONFIG, NOISE_SCALARS, and
# get_calibration_prior(). Run with:
#   CALIBRATION_CONFIG=lai.jl bash experiments/calibration/run_calibration.sh
#
# This calibrates the prognostic optimal-LAI model (Zhou et al. 2025,
# `ZhouOptimalLAIModel`) so the forward model is run with `prognostic_lai = true`
# (canopy computes LAI from P-model potential GPP instead of prescribing MODIS).
# The single target is the model `lai` diagnostic (m^2 m^-2) compared against the
# MODIS LAI observation (`get_modis_lai_obs_var` in data_sources.jl), as seasonal
# averages.
#
# Ensemble size = 4 params * 2 + 1 = 9 (TransformUnscented); set the number of
# tasks in run_calibration.sh / the batch script to 9.

"""Per-variable observation variance; multiplies the identity in `ScalarCovariance`.
Units are (m^2 m^-2)^2.

A scalar of 0.5 corresponds to an observation σ ≈ 0.7 m^2 m^-2, which is a
reasonable combined MODIS-retrieval + model-structural error for seasonal LAI
(typical LAI ranges 0–6). Loosen it (→ 1.0+) if EKI overfits the seasonal
amplitude at the expense of timing; tighten it (→ 0.25) to pull harder on the
peak-LAI magnitude.
"""
const NOISE_SCALARS = Dict("lai" => 0.5)

"""
    get_calibration_prior()

Prior over the 4 optimal-LAI parameters that shape the simulated LAI seasonal
cycle and magnitude. Means are the TOML defaults; `optimal_lai_k` (light
extinction, 0.5) is left fixed because it is strongly correlated with
`optimal_lai_z` (both enter the energy-limited fAPAR = 1 - z/(k*A0)) and would
otherwise make the inversion degenerate.

TransformUnscented ensemble = 4*2+1 = 9.

  - optimal_lai_z     leaf construction/maintenance cost (mol m^-2 yr^-1). Sets
                      the energy-limited LAI_max via 1 - z/(k*A0): larger z ->
                      lower peak LAI.
  - optimal_lai_sigma departure from square-wave LAI dynamics (dimensionless).
                      Controls the m ratio (Eq. 20) and hence how sharply LAI
                      tracks the growing season.
  - optimal_lai_alpha EMA smoothing factor (dimensionless). Sets the
                      green-up/senescence lag (alpha ~ 0.067 ≈ 15-day memory).
  - optimal_lai_f0    fraction of precipitation available for transpiration
                      (dimensionless). Sets the water-limited LAI_max in
                      moisture-limited regions.
"""
function get_calibration_prior()
    priors = [
        EKP.constrained_gaussian("optimal_lai_z", 12.227, 4.0, 1.0, 40.0),
        EKP.constrained_gaussian("optimal_lai_sigma", 1.1, 0.3, 0.1, 3.0),
        # alpha and f0 are bounded fractions; σ kept clear of the bounds.
        EKP.constrained_gaussian("optimal_lai_alpha", 0.067, 0.03, 0.01, 0.3),
        EKP.constrained_gaussian("optimal_lai_f0", 0.65, 0.2, 0.05, 1.0),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lai"],
    minibatch_size = 1,
    n_iterations = 10,
    # One annual cycle to calibrate on (DJF-MAM-JJA from Dec 1 -> Sep 1, plus
    # SON via extend = Month(3)). With spinup = Year(2), the simulation runs
    # 2015-12-01 -> 2018-12-01: two years of spin-up so the optimal-LAI annual
    # potential GPP (A0_annual) self-corrects off its stale file value, then one
    # year scored against MODIS. All dates stay within MODIS (2000-2020) and the
    # high-res ERA5 forcing coverage. Add more (year) tuples here (and raise
    # minibatch_size) to calibrate across multiple years.
    sample_date_ranges = [("2017-12-1", "2018-9-1")],
    extend = Dates.Month(3),
    spinup = Dates.Year(2),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_lai.jld2",
    model_type = ClimaLand.LandModel,
    prognostic_lai = true,
)
