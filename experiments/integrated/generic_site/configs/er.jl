# Calibration of soil heterotrophic respiration (DAMM kinetics) at a single
# fluxnet site. Targets ecosystem respiration (`er`) compared against the
# fluxnet `RECO_NT_VUT_REF` column.
#
# Loaded by experiments/integrated/generic_site/calibrate_site.jl. Defines:
#   CALIBRATE_CONFIG, NOISE_VARIANCES, get_calibration_prior()
#
# Loss aggregation: daily means over each yearly sample window. Each sample
# is a 365-day window preceded by a 90-day spinup. Up to `max_n_samples`
# windows per site (sites with ≥5 years are capped; sites with fewer fall
# back to whatever they have). Minibatching draws `minibatch_size` samples
# per EKP iteration; with `minibatch_size = 1`, each iteration runs one year.

const CALIBRATE_CONFIG = (;
    short_names = ["er"],
    n_iterations = 5,
    cal_window_days = 365,
    spinup_days = 90,
    minibatch_size = 1,
    max_n_samples = 5,
    rng_seed = 42,
)

# (g C m^-2 day^-1)^2 per daily bin. σ ≈ 2 g C m^-2 day^-1 — fluxnet RECO
# typical range ~0–10 g C m^-2 day^-1, so this is ~20% relative.
const NOISE_VARIANCES = Dict("er" => 4.0)

"""
    get_calibration_prior()

Four DAMM kinetic parameters mirroring the DAMM block of
`experiments/calibration/configs/gpp_er_nee_lhf.jl`, with
`O2_michaelis_constant` widened so μ−σ stays well above zero (the reference
prior leaves μ−σ ≈ 5e-4, near the lower bound, which makes EKI lean against
the bound and converge poorly on single-site data).
"""
function get_calibration_prior()
    return EKP.combine_distributions([
        EKP.constrained_gaussian(
            "soilCO2_reference_rate",
            2.17564e-7,
            2.0e-7,
            1.0e-9,
            2.0e-6,
        ),
        EKP.constrained_gaussian(
            "soilCO2_activation_energy",
            37357.0,
            25000.0,
            5000.0,
            150000.0,
        ),
        EKP.constrained_gaussian(
            "michaelis_constant",
            0.459997,
            0.25,
            0.001,
            2.0,
        ),
        # Widened relative to gpp_er_nee_lhf.jl (μ 0.00255 → 0.05, σ 0.002 → 0.02).
        EKP.constrained_gaussian(
            "O2_michaelis_constant",
            0.05,
            0.02,
            0.005,
            0.2,
        ),
    ])
end
