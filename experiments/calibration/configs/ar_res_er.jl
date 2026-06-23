# Calibration config: autotrophic respiration vs inversion residual ER.
# Included by run_calibration.jl; defines CALIBRATE_CONFIG, NOISE_SCALARS, and
# get_calibration_prior(). Run with:
#   CALIBRATION_CONFIG=ar_res_er.jl \
#     bash experiments/calibration/run_calibration.sh
#
# Goal: pull total ecosystem respiration ER = Ra + Rh onto the inversion
# residual ER (res_er = CT2022 NEE + GOSIF GPP) by moving ONLY autotrophic
# respiration. GPP and the DAMM/Rh params are already baked into
# toml/default_parameters.toml (GOSIF/energy + Hashimoto calibrations), so Rh is
# held fixed here and ER is shifted entirely through Ra.
#
# The two free parameters live solely in autotrophic_respiration.jl and never
# enter the photosynthesis / stomatal-conductance / GPP path, so GPP is
# unchanged by this calibration (autotrophic respiration is computed strictly
# downstream of photosynthesis: it reads An_canopy/Rd_canopy and returns Ra,
# and nothing reads Ra back; GPP is reported as gross assimilation):
#   - autotrophic_respiration_Rd_ref      base root/stem maintenance
#                                          (R_root = Rd_ref*RAI, R_stem = Rd_ref*μs*SAI)
#   - relative_contribution_factor (Rel)  growth respiration (Rg = Rel*max(An-Rpm,0))
# The leaf maintenance term R_leaf = Rd_canopy is inherited from the (fixed)
# photosynthesis Rd, so it is not a lever here.
#
# Target (from the `inversion_nee` artifact, g C m^-2 day^-1, + = source):
#   res_er  ER residual = NEE + GPP. observation_map.jl reads the model `er`
#           diagnostic and retags it to res_er.

"""Per-variable observation variance; multiplies the identity in `ScalarCovariance`.
Units are (dataset units)^2: res_er in (g C m^-2 day^-1)^2.

Matches the res_er weight used in inversionnee_sifgpp_er_lhf_shf_lwu.jl; res_er
is a residual (NEE + GPP) so it carries real structural noise.
"""
const NOISE_SCALARS = Dict("res_er" => 0.5)

"""
    get_calibration_prior()

Prior over the two autotrophic-respiration parameters;
TransformUnscented ensemble = 2*2+1 = 5 members.

Both priors are centered on the current `toml/default_parameters.toml` defaults.
autotrophic_respiration_Rd_ref uses a [0, Inf] lognormal transform (a finite
upper bound makes constrained_gaussian silently fall back to the bound midpoint
for these small-magnitude rates; see calibration_soilco2_prior_misspecified) and
is given a broad σ so EKI can lower Ra (the model over-respires) without hitting
a bound. relative_contribution_factor is two-sided on [0, 1].
"""
function get_calibration_prior()
    priors = [
        EKP.constrained_gaussian(
            "autotrophic_respiration_Rd_ref",
            7.83e-7,
            5.0e-7,
            0.0,
            Inf,
        ),
        EKP.constrained_gaussian(
            "relative_contribution_factor",
            0.3,
            0.15,
            0.0,
            1.0,
        ),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["res_er"],
    minibatch_size = 2,
    n_iterations = 10,
    # Yearly DJF-MAM-JJA samples (Dec 1 -> Sep 1). With spinup = Year(1) and
    # extend = Month(3), all 16 samples stay within the artifact's 2002-2020
    # coverage; 16 is divisible by minibatch_size so none is dropped per epoch.
    sample_date_ranges = [
        ("$(year)-12-1", "$(year+1)-9-1") for year in 2003:2018
    ],
    extend = Dates.Month(3),
    # 1-year spinup so deep soil temperature/moisture settle before ER vs obs.
    spinup = Dates.Year(1),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_ar_res_er.jld2",
    model_type = ClimaLand.LandModel,
    use_rosetta = true,
)
