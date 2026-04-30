# Default calibration configuration for single-site fluxnet calibration.
#
# Loaded by experiments/integrated/generic_site/calibrate_site.jl, which
# expects the following names to be defined:
#
#   - get_calibration_prior() -> EKP combined distribution
#   - CALIBRATE_CONFIG        -> NamedTuple with keys
#         short_names         :: Vector{String}   variables in the obs vector
#         n_iterations        :: Int              EKP iterations
#         cal_duration_days   :: Int              calibration window length
#         rng_seed            :: Int
#   - NOISE_VARIANCES         -> Dict mapping each short_name in
#                                 CALIBRATE_CONFIG.short_names to a per-bin
#                                 scalar variance used to build the diagonal
#                                 noise covariance.
#
# To create a new config, copy this file (e.g. configs/my_experiment.jl) and
# adjust the settings below. Select it at run time with the env var
# CALIBRATION_CONFIG, e.g. CALIBRATION_CONFIG=my_experiment.jl julia ...

const CALIBRATE_CONFIG = (;
    short_names = ["gpp", "lhf"],
    n_iterations = 5,
    cal_duration_days = 365,
    rng_seed = 42,
)

# Per-variable scalar variance applied to every monthly bin. GPP scale
# ~3 µmol m⁻² s⁻¹ = 3e-6 mol m⁻² s⁻¹; LHF scale ~30 W m⁻².
const NOISE_VARIANCES = Dict("gpp" => (3e-6)^2, "lhf" => (30.0)^2)

"""
    get_calibration_prior()

Two PModel-related parameters, mirroring a subset of
`experiments/calibration/configs/gpp_lhf.jl`. Easy to expand as the calibration
pipeline matures.
"""
function get_calibration_prior()
    return EKP.combine_distributions([
        EKP.constrained_gaussian("pmodel_cstar", 0.41, 0.05, 0.2, 0.7),
        EKP.constrained_gaussian("moisture_stress_c", 0.27, 0.15, 0.05, 1.0),
    ])
end
