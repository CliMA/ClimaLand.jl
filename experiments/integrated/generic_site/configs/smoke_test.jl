# Short-window calibration config for CI smoke tests.
#
# Same prior and short_names as configs/default.jl but with a 90-day window
# and 2 iterations so a full ensemble run finishes inside a CI step.

const CALIBRATE_CONFIG = (;
    short_names = ["gpp", "lhf"],
    n_iterations = 2,
    cal_window_days = 90,
    spinup_days = 0,
    minibatch_size = 1,
    max_n_samples = 1,
    rng_seed = 42,
)

const NOISE_VARIANCES = Dict("gpp" => (3e-6)^2, "lhf" => (30.0)^2)

function get_calibration_prior()
    return EKP.combine_distributions([
        EKP.constrained_gaussian("pmodel_cstar", 0.41, 0.05, 0.2, 0.7),
        EKP.constrained_gaussian("moisture_stress_c", 0.27, 0.15, 0.05, 1.0),
    ])
end
