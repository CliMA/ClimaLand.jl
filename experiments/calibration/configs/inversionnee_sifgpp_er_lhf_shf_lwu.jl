# Calibration config: inversion NEE/GPP/ER/Rh + ERA5 LHF/SHF/LWU.
# Included by run_calibration.jl; defines CALIBRATE_CONFIG, NOISE_SCALARS, and
# get_calibration_prior(). Run with:
#   CALIBRATION_CONFIG=inversionnee_sifgpp_er_lhf_shf_lwu.jl \
#     bash experiments/calibration/run_calibration.sh
#
# Targets use distinct short_names so they don't collide with the FLUXCOM
# gpp/er/nee keys. Carbon targets come from the `inversion_nee` artifact in
# g C m^-2 day^-1 (sign: + = source, except GPP + = uptake):
#   inv_nee  NEE, CarbonTracker CT2022
#   sif_gpp  GPP, GOSIF-GPP v2
#   res_er   ER residual = NEE + GPP
#   inv_hr   Rh, Hashimoto 2015 (2013-2020 filled with 2002-2012 climatology)
# plus ERA5 lhf/shf/lwu in W m^-2 (+ = upward). Sign conventions match the
# model, which computes NEE = ER - GPP. data_sources.jl converts nee/gpp/er
# from monthly totals to daily rates using actual days-in-month (not
# 365.25/12); Rh is already daily and not rescaled.

"""Per-variable observation variances; multiply the identity in `ScalarCovariance`.
Units are (dataset units)^2: carbon in (g C m^-2 day^-1)^2, energy in (W m^-2)^2.

inv_hr is the tightest carbon target (0.05) so soilCO2 params are pulled to the
Hashimoto magnitude instead of collapsing. inv_nee is loose (0.2) so EKI fits
NEE's spatial structure without buying amplitude via an unphysical AR/Rh split
(the two errors cancel in NEE). res_er ≡ inv_nee + sif_gpp, so it carries no
independent info; kept moderate (0.5). GPP stays loose (2.0) to absorb the
P-model high-latitude bias. Energy variances sit at the RMSE floor (LWU 25,
SHF/LHF 144) and act as background constraints.
"""
const NOISE_SCALARS = Dict(
    "inv_nee" => 0.2,
    "sif_gpp" => 2.0,
    "res_er" => 0.5,
    "inv_hr" => 0.05,
    "lwu" => 25.0,
    "shf" => 144.0,
    "lhf" => 144.0,
)

"""
    get_calibration_prior()

Combined prior over 15 parameters (5 P-model/moisture-stress + 4 DAMM Rh + 1
JULES AR Rd_ref + 5 canopy transfer); TransformUnscented ensemble = 15*2+1 = 31.

relative_contribution_factor, autotrophic_respiration_Q10 and
root_leaf_nitrogen_ratio (μr) are no longer calibrated (fixed in TOML): they
were the AR-inflation / partition-degeneracy levers — with them free, EKI bought
NEE amplitude by inflating AR and collapsing Rh. Means are recentered on the
eki_file.jld2 iter-3 ensemble mean, except the canopy-transfer params (whose
iter trajectory thrashed) which keep physical defaults. soilCO2_reference_rate,
O2_michaelis_constant and Rd_ref use a [0, Inf] lognormal transform; finite-bound
constrained_gaussian silently fell back to the bound midpoint for these.
"""
function get_calibration_prior()
    priors = [
        # P-model + moisture stress (means = eki_file iter-3 ensemble mean)
        EKP.constrained_gaussian("pmodel_cstar", 0.4335786422686187, 0.05, 0.2, 0.7),
        EKP.constrained_gaussian("pmodel_β_c3", 85.06089150945417, 40.0, 10.0, 300.0),
        EKP.constrained_gaussian("pmodel_β_c4", 28.65798918466436, 12.0, 5.0, 100.0),
        # σ kept tight so μ+σ stays clear of the upper bound 0.999.
        EKP.constrained_gaussian("pmodel_α", 0.9543827460568335, 0.012, 0.85, 0.999),
        EKP.constrained_gaussian(
            "moisture_stress_c",
            0.505925674083973,
            0.15,
            0.05,
            1.0,
        ),
        # DAMM heterotrophic respiration.
        # soilCO2_reference_rate ([0, Inf] lognormal): 6e-8 gives 4-site summer
        # Rh 1.84 vs obs 1.80 in the single-column sweep (priormean_sweep.jl).
        # Rh depends only on this rate; activation_energy refines its spatial spread.
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
        # autotrophic_respiration_Rd_ref ([0, Inf] lognormal): base rate for the
        # LAI-independent root/stem maintenance (R_root = Rd_ref*μr*RAI*β,
        # R_stem = Rd_ref*μs*(ηsl*h/σl)*SAI). 2e-7 gives the best summer AR in the
        # single-column sweep with μr fixed at 1.0 (grassland 2.86 vs obs 2.74).
        EKP.constrained_gaussian(
            "autotrophic_respiration_Rd_ref",
            2.0e-7,
            1.5e-7,
            0.0,
            Inf,
        ),
        # Canopy turbulent / radiative transfer. Kept at physical means (iter-3
        # is a poor center: e.g. canopy_d_coeff thrashed to ~0.01).
        EKP.constrained_gaussian("leaf_Cd", 0.08, 0.03, 0, Inf),
        EKP.constrained_gaussian("canopy_z_0m_coeff", 0.13, 0.05, 0.0, 0.5),
        EKP.constrained_gaussian("canopy_z_0b_coeff", 0.013, 0.005, 0.0, 0.05),
        EKP.constrained_gaussian("canopy_d_coeff", 0.67, 0.2, 0.0, 1.0),
        EKP.constrained_gaussian("canopy_K_lw", 1.0, 0.2, 0.0, 2.0),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["inv_nee", "sif_gpp", "res_er", "inv_hr", "lhf", "shf", "lwu"],
    minibatch_size = 2,
    n_iterations = 10,
    # Yearly DJF-MAM-JJA samples (Dec 1 -> Sep 1). With spinup = Year(1) and
    # extend = Month(3), all 16 samples stay within the artifact's 2002-2020
    # coverage; 16 is divisible by minibatch_size so none is dropped per epoch.
    sample_date_ranges = [
        ("$(year)-12-1", "$(year+1)-9-1") for year in 2003:2018
    ],
    extend = Dates.Month(3),
    # 1-year spinup so deep soil temperature/moisture settle before Hr vs ER.
    spinup = Dates.Year(1),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_inversionnee_sifgpp_er_lhf_shf_lwu.jld2",
    model_type = ClimaLand.LandModel,
)
