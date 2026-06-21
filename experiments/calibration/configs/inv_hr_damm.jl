# Calibration config: DAMM heterotrophic respiration only.
# Included by run_calibration.jl; defines CALIBRATE_CONFIG, NOISE_SCALARS, and
# get_calibration_prior(). Run with:
#   CALIBRATION_CONFIG=inv_hr_damm.jl \
#     bash experiments/calibration/run_calibration.sh
#
# This isolates the soil-CO2 (DAMM) heterotrophic-respiration parameters: it
# calibrates the four DAMM parameters against a single Rh target and nothing
# else. The autotrophic-respiration and P-model/canopy parameters that the
# larger inversionnee_sifgpp_er_lhf_shf_lwu.jl config also moves are excluded,
# since Rh alone cannot constrain them.
#
# Target (from the `inversion_nee` artifact, g C m^-2 day^-1, + = source):
#   inv_hr  Rh, Hashimoto 2015 (2013-2020 filled with 2002-2012 climatology).
# observation_map.jl reads the model `hr` diagnostic and retags it to inv_hr.
# Rh is already a daily rate, so it is not rescaled from monthly totals.

"""Per-variable observation variance; multiplies the identity in `ScalarCovariance`.
Units are (dataset units)^2: Rh in (g C m^-2 day^-1)^2.

inv_hr is kept tight (0.05, matching inversionnee_sifgpp_er_lhf_shf_lwu.jl) so
the DAMM rate is pulled to the Hashimoto magnitude rather than drifting.
"""
const NOISE_SCALARS = Dict("inv_hr" => 0.05)

"""
    get_calibration_prior()

Prior over the four DAMM heterotrophic-respiration parameters;
TransformUnscented ensemble = 4*2+1 = 9 members.

Means, σ and bounds are inherited from the DAMM block of
inversionnee_sifgpp_er_lhf_shf_lwu.jl. soilCO2_reference_rate and
O2_michaelis_constant use a [0, Inf] lognormal transform: a finite upper bound
makes constrained_gaussian silently fall back to the bound midpoint for these
small-magnitude rates (see calibration_soilco2_prior_misspecified). Rh depends
chiefly on soilCO2_reference_rate; soilCO2_activation_energy refines the spatial
spread, and the two Michaelis constants shape the moisture/oxygen limitation.
"""
function get_calibration_prior()
    priors = [
        # soilCO2_reference_rate ([0, Inf] lognormal): 6e-8 gives 4-site summer
        # Rh 1.84 vs obs 1.80 in the single-column sweep (priormean_sweep.jl).
        EKP.constrained_gaussian(
            "soilCO2_reference_rate",
            6.0e-8,
            4.0e-8,
            0.0,
            Inf,
        ),
        EKP.constrained_gaussian(
            "soilCO2_activation_energy",
            36897.373510987745,
            20000.0,
            5000.0,
            150000.0,
        ),
        # σ reduced to 0.03 so μ-σ stays above the lower bound 0.001.
        EKP.constrained_gaussian(
            "michaelis_constant",
            0.05113112510036322,
            0.03,
            0.001,
            2.0,
        ),
        # [0, Inf] lognormal (same fallback issue as soilCO2_reference_rate).
        EKP.constrained_gaussian(
            "O2_michaelis_constant",
            1.9e-4,
            1.5e-4,
            0.0,
            Inf,
        ),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["inv_hr"],
    minibatch_size = 2,
    n_iterations = 10,
    # Yearly DJF-MAM-JJA samples (Dec 1 -> Sep 1). With spinup = Year(1) and
    # extend = Month(3), all 16 samples stay within the artifact's 2002-2020
    # coverage; 16 is divisible by minibatch_size so none is dropped per epoch.
    sample_date_ranges = [
        ("$(year)-12-1", "$(year+1)-9-1") for year in 2003:2018
    ],
    extend = Dates.Month(3),
    # 1-year spinup so deep soil temperature/moisture settle before Hr vs obs.
    spinup = Dates.Year(1),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_inv_hr_damm.jld2",
    model_type = ClimaLand.LandModel,
    use_rosetta = true,
)
