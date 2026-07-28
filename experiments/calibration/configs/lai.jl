# Calibration config: optimal-LAI parameters against MODIS LAI.
# Included by run_calibration.jl; defines CALIBRATE_CONFIG, NOISE_SCALARS, and
# get_calibration_prior(). Run with:
#   CALIBRATION_CONFIG=lai.jl bash experiments/calibration/run_calibration.sh
#
# This calibrates the prognostic optimal-LAI model (Zhou et al. 2025,
# `ZhouOptimalLAIModel`), so the forward model runs with `prognostic_lai = true`
# and the canopy computes LAI from the P-model potential GPP instead of
# prescribing MODIS. The single target is the model `lai` diagnostic (m^2 m^-2)
# against the MODIS LAI observation (`get_modis_lai_obs_var` in data_sources.jl),
# as seasonal averages.
#
# Stage 3 of the LAI calibration pipeline. It depends on both earlier stages:
#   1. the P-model/energy calibration (sifgpp_lhf_shf_lwu_rosetta.jl) sets the
#      potential GPP that drives LAI, so run it first and adopt its parameters;
#   2. the rebuilt optimal-LAI initial-condition artifact means the trailing
#      totals start near equilibrium, which is what lets `spinup` stay short
#      (see the note on spinup below).
#
# Ensemble size = 5 params * 2 + 1 = 11 (TransformUnscented); set the number of
# tasks in run_calibration.sh / the batch script to 11.

"""Per-variable observation variance; multiplies the identity in `ScalarCovariance`.
Units are (m^2 m^-2)^2.

0.5 corresponds to an observation σ ≈ 0.7 m^2 m^-2, a reasonable combined
MODIS-retrieval + model-structural error for seasonal LAI (typical range 0-6).

This is a tuning knob, not a measured quantity. Loosen it (-> 1.0+) if EKI
overfits the seasonal amplitude at the expense of timing; tighten it (-> 0.25) to
pull harder on peak-LAI magnitude. With the C3/C4 split below there are five
parameters rather than three, so there is more freedom to overfit; if the
posterior collapses onto extreme z/sigma values, loosen this before widening the
priors.
"""
const NOISE_SCALARS = Dict("lai" => 0.5)

"""
    get_calibration_prior()

Prior over the 5 optimal-LAI parameters that shape the simulated LAI seasonal
cycle and magnitude. Means are the current TOML defaults.

`optimal_lai_k` (light extinction, 0.5) is left fixed: it is strongly correlated
with `optimal_lai_z` (both enter the energy-limited fAPAR = 1 - z/(k*A0)) and
would make the inversion degenerate. `optimal_lai_f0` is excluded because it is
inert — f0 is computed each step from the aridity index PET_annual/precip_annual,
so the TOML scalar is never read. `optimal_lai_z_a0` is left at 0 (uniform z in
A0) to keep this inversion identifiable; calibrate it separately if the wet
tropics stay biased after this stage.

TransformUnscented ensemble = 5*2+1 = 11.

  - optimal_lai_z / z_c4         leaf construction+maintenance cost
                                 (mol m^-2 yr^-1), C3 and C4. Sets the
                                 energy-limited LAI_max via 1 - z/(k*A0):
                                 larger z -> lower peak LAI. Blended per column
                                 by the model's own C3 fraction (`pft_blend`),
                                 so the split only bites where C4 is present.
  - optimal_lai_sigma / sigma_c4 departure from square-wave LAI dynamics
                                 (dimensionless), C3 and C4. Controls how
                                 sharply LAI tracks the growing season.
  - optimal_lai_alpha            LAI acclimation rate (dimensionless); the
                                 green-up/senescence lag is τ = 1 day / alpha,
                                 so 0.0701 ≈ 14 days. Not split by pathway: it
                                 enters as the `RunningMean` timescale of a
                                 `TimeIntegratedVariable`, which is a scalar for
                                 the whole field and cannot vary per column.

The C3 and C4 priors are deliberately identical: the defaults are equal, so the
prior expresses no belief about which pathway is costlier, and any split the
posterior finds is driven by the MODIS data.
"""
function get_calibration_prior()
    priors = [
        EKP.constrained_gaussian("optimal_lai_z", 15.0, 4.0, 1.0, 40.0),
        EKP.constrained_gaussian("optimal_lai_z_c4", 15.0, 4.0, 1.0, 40.0),
        EKP.constrained_gaussian("optimal_lai_sigma", 0.939, 0.3, 0.1, 3.0),
        EKP.constrained_gaussian("optimal_lai_sigma_c4", 0.939, 0.3, 0.1, 3.0),
        # alpha is a bounded rate; σ kept clear of the bounds.
        EKP.constrained_gaussian("optimal_lai_alpha", 0.0701, 0.03, 0.01, 0.3),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lai"],
    minibatch_size = 1,
    n_iterations = 10,
    # One annual cycle to score (DJF-MAM-JJA from Dec 1 -> Sep 1, plus SON via
    # extend = Month(3)). Add more (year) tuples here, and raise minibatch_size,
    # to calibrate across multiple years.
    sample_date_ranges = [("2017-12-1", "2018-9-1")],
    extend = Dates.Month(3),
    # `optimal_lai_tau_long_term` is 2 years, so the trailing totals need roughly
    # that long to forget a bad initial condition. One year of spinup is only
    # enough because stage 2 rebuilds the IC artifact from an equilibrated
    # 10-year run — with the stock artifact this must go back up to Year(3) or
    # the inversion fits spin-up transients instead of the seasonal cycle.
    spinup = Dates.Year(1),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_lai.jld2",
    model_type = ClimaLand.LandModel,
    prognostic_lai = true,
)
