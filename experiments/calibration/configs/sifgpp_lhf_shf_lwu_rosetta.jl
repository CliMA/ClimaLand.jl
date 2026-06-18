# Calibration config: GOSIF GPP + ERA5 LHF/SHF/LWU, with a Rosetta soil.
# Included by run_calibration.jl; defines CALIBRATE_CONFIG, NOISE_SCALARS, and
# get_calibration_prior(). Run with:
#   CALIBRATION_CONFIG=sifgpp_lhf_shf_lwu_rosetta.jl \
#     bash experiments/calibration/run_calibration.sh
#
# This is the energy/carbon-uptake subset of inversionnee_sifgpp_er_lhf_shf_lwu.jl:
# the carbon-source targets (inv_nee/res_er/inv_hr) are dropped, so only the
# parameters that GPP/LHF/SHF/LWU can actually constrain are calibrated. The
# soilCO2/DAMM-Rh params and autotrophic_respiration_Rd_ref are therefore not
# included (they only move Rh/AR, which have no target here).
#
# Targets:
#   sif_gpp  GPP, GOSIF-GPP v2 (from the `inversion_nee` artifact, + = uptake)
#   lhf/shf/lwu  ERA5 turbulent + upwelling longwave fluxes in W m^-2 (+ = upward)
#
# The model uses the Rosetta (Montzka et al. 2017) van Genuchten soil retention
# parameters via `use_rosetta = true` (the default), which also selects the
# matching spun-up Rosetta initial conditions (`rosetta_spunup_ic`, PR #1772) so
# the soil state is consistent with the retention curve.

"""Per-variable observation variances; multiply the identity in `ScalarCovariance`.
Units are (dataset units)^2: GPP in (g C m^-2 day^-1)^2, energy in (W m^-2)^2.

GPP is tightened to 0.6 (vs the inversion config's 2.0) since it is the only
carbon target and we want the P-model params pulled to GOSIF. Energy variances
sit near the RMSE floor (LWU 25, SHF/LHF 100) and act as background constraints.
"""
const NOISE_SCALARS = Dict(
    "sif_gpp" => 0.6,
    "lhf" => 100.0,
    "shf" => 100.0,
    "lwu" => 25.0,
)

"""
    get_calibration_prior()

Prior over the 10 parameters relevant to GPP/LHF/SHF/LWU (5 P-model/moisture
stress + 5 canopy turbulent/radiative transfer); TransformUnscented ensemble =
10*2+1 = 21.

Each prior mean equals the corresponding value in `toml/default_parameters.toml`
(synced to the `ar/calibrated_C_energy_together` branch), so the prior is
centered on the model defaults. The σ and bounds are inherited from the relevant
subset of `inversionnee_sifgpp_er_lhf_shf_lwu.jl`, except where the calibrated
default sat at a former bound — there the bound is relaxed and σ shrunk so the
constrained_gaussian can actually center on the default (EKP errors when
μ ± σ crosses a bound):
  - pmodel_β_c4:      σ 12 -> 2,    lower 5 -> 1     (default 5.36 was at lower=5)
  - canopy_z_0b_coeff: σ 0.005 -> 0.01, upper 0.05 -> 0.1 (default 0.044 near upper=0.05)
  - canopy_d_coeff:    σ 0.2 -> 0.025                 (default 0.057 near lower=0)
"""
function get_calibration_prior()
    priors = [
        # P-model + moisture stress. Constrain GPP and, via stomatal
        # conductance, transpiration (LHF).
        EKP.constrained_gaussian("pmodel_cstar", 0.42890384025557793, 0.05, 0.2, 0.7),
        EKP.constrained_gaussian("pmodel_β_c3", 91.09351405009605, 40.0, 10.0, 300.0),
        EKP.constrained_gaussian("pmodel_β_c4", 5.356145477815346, 2.0, 1.0, 100.0),
        # σ kept tight so μ+σ stays clear of the upper bound 0.999.
        EKP.constrained_gaussian("pmodel_α", 0.9723018635387689, 0.012, 0.85, 0.999),
        EKP.constrained_gaussian(
            "moisture_stress_c",
            0.591401756353606,
            0.15,
            0.05,
            1.0,
        ),
        # Canopy turbulent / radiative transfer. These constrain SHF/LHF
        # (roughness/drag) and LWU (longwave extinction).
        EKP.constrained_gaussian("leaf_Cd", 0.07257315811847019, 0.03, 0, Inf),
        EKP.constrained_gaussian("canopy_z_0m_coeff", 0.3494203357014585, 0.05, 0.0, 0.5),
        EKP.constrained_gaussian("canopy_z_0b_coeff", 0.04435103776675633, 0.01, 0.0, 0.1),
        EKP.constrained_gaussian("canopy_d_coeff", 0.05725798250714986, 0.025, 0.0, 1.0),
        EKP.constrained_gaussian("canopy_K_lw", 0.9191833718378635, 0.2, 0.0, 2.0),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["sif_gpp", "lhf", "shf", "lwu"],
    minibatch_size = 2,
    n_iterations = 10,
    # Yearly DJF-MAM-JJA samples (Dec 1 -> Sep 1). With spinup = Year(1) and
    # extend = Month(3), all 16 samples stay within the artifact's 2002-2020
    # coverage; 16 is divisible by minibatch_size so none is dropped per epoch.
    sample_date_ranges = [
        ("$(year)-12-1", "$(year+1)-9-1") for year in 2003:2018
    ],
    extend = Dates.Month(3),
    # 1-year spinup so soil moisture/temperature settle before the flux average.
    spinup = Dates.Year(1),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_sifgpp_lhf_shf_lwu.jld2",
    model_type = ClimaLand.LandModel,
    use_rosetta = true,
)
