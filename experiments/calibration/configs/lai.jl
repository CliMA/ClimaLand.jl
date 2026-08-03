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
# Ensemble size = 5 params * 2 + 1 = 11 (TransformUnscented). On the Derecho
# backend each member gets its own single-task PBS job, so nothing needs setting
# by hand; expect 11 member jobs per iteration.
#
# The `lai` target is masked to natural, undisturbed vegetation — see
# build_natural_vegetation_mask.jl and apply_natural_vegetation_mask in
# generate_observations.jl. Cropland, urban, glacier and habitually-burning
# cells are excluded because this model represents none of those processes.

"""Per-variable observation variance; multiplies the identity in `ScalarCovariance`.
Units are (m^2 m^-2)^2.

0.5 corresponds to an observation σ ≈ 0.7 m^2 m^-2, a reasonable combined
MODIS-retrieval + model-structural error for seasonal LAI (typical range 0-6).

⚠ MEASURED TO BE A NO-OP FOR THIS CONFIG. Halving it to 0.25 scaled the misfit
by exactly the variance ratio at every iteration (2.00, 2.01, 2.01, 1.99) while
the parameters tracked the 0.5 run to within 1-2 %, and the final masked LAI RMSE
was 0.883 against 0.885. `DataMisfitController` sets its timestep from the misfit
magnitude, so uniformly rescaling the covariance is divided straight back out.

The scalar DOES matter for the RELATIVE weighting between several targets — that
is how `sifgpp_lhf_shf_lwu_rosetta.jl` uses it across four variables. It does not
give a lever on amplitude-versus-timing in a SINGLE-target calibration under an
adaptive-timestep scheduler. Do not reach for it here expecting one; the binding
constraint is the z-sigma degeneracy, not the assumed observation error.
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
        # alpha is a bounded rate; σ kept clear of the bounds. Do not narrow the
        # upper bound to dodge the GPU DomainError that kills ~2.5 % of members:
        # over 275 member-runs alpha does not separate the failures (crashes at
        # 0.122, successes at 0.234), and since every calibration settles at
        # 0.18-0.21 a cap there binds on the posterior. See the census and the
        # crash diagnosis in unsupervised_loop/STATE.md.
        EKP.constrained_gaussian("optimal_lai_alpha", 0.0701, 0.03, 0.01, 0.3),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lai"],
    minibatch_size = 1,
    # N_ITERATIONS lets a run be advanced one iteration at a time. ClimaCalibrate
    # resumes at `last_completed_iteration + 1` and stops at `n_iterations`, so
    # raising it in steps keeps each orchestrator job short enough to exit
    # normally inside a queue's walltime cap, instead of being killed mid-
    # iteration and leaving that iteration's member jobs running as orphans that
    # race the resubmitted orchestrator for the same member directories.
    # Default is unchanged.
    n_iterations = parse(Int, get(ENV, "N_ITERATIONS", "5")),
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
