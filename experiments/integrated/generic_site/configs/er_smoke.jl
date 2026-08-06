# CI smoke variant of `er.jl`: shorter window and fewer iterations so the
# full ensemble run finishes inside a buildkite step.

const CALIBRATE_CONFIG = (;
    short_names = ["er"],
    n_iterations = 1,
    cal_window_days = 30,
    spinup_days = 15,
    minibatch_size = 1,
    max_n_samples = 1,
    rng_seed = 42,
)

const NOISE_VARIANCES = Dict("er" => 4.0)

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
        EKP.constrained_gaussian(
            "O2_michaelis_constant",
            0.05,
            0.02,
            0.005,
            0.2,
        ),
    ])
end
