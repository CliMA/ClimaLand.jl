# TODO: I am still iffy about this function, because the only reason why I need
# it is to be able to use it in the worker processes. For example, I use
# sample_date_ranges, but for some fields like minibatch_size and n_iterations,
# I only need to know them at the start to initialize the variables. Another
# thing that I can do is to divide them into two or more functions, by what is
# used in the workers and observation making and what is needed in just the
# iterations

# TODO: Add accessor functions in ClimaAnalysis, so that something can be
# implemented in ClimaCalibrate to check for units, dates, and grid

# TODO: Don't know if importing packages in workers cause increase in compile times

# Intead of using constants, we use functions to be able to pass arguments to
# the workers.

import Dates
import Random
import EnsembleKalmanProcesses as EKP
import JLD2
import ClimaCalibrate

"""
    get_config()

Return a named tuple consisting of
- `nelements`: Resolution of the model. Used to determine the correct resolution
  for the mask.
- `short_names`: Short names of the variables that are included in the
  calibration.
- `sample_date_ranges`: Tuples of start and end dates for each of the samples.
- `spinup`: The amount to spin up the simulation for.
- `extend`: The amount to extend the simulation for from the last date of the
  minibatch.
- `output_dir`: The directory to store the iterations and members of the
  calibration.

!!! note "Observations"
    If `short_names` or `sample_date_ranges` change, then you need to remake the
    observations
"""
function get_config()
    return (;
        nelements = (101, 15),
        short_names = ["swu"], # ["lhf", "shf", "swu", "lwu"],

        # Time ranges of each sample which is typically is a single year. Because we are
        # dealing with seasons, we use December and September because those are the
        # start of DJF and SON
        sample_date_ranges = [
            (Dates.DateTime(year, 12, 1), Dates.DateTime(year + 1, 9, 1)) for
            year in 2000:2009
        ],
        spinup = Dates.Month(3), # Dates.Year(1)
        # TODO: Not sure if sample_date_ranges is the best way of determining
        # the samples

        # The `extend = Dates.Month(3)` comes from calibrating on seasonal averages
        # We want a full season and since the date of a season should be the first
        # date of the season, we want to do 3 more months to capture the full season
        extend = Dates.Month(3),
        output_dir = "experiments/better_calibration/land_model",
    )
end

"""
    get_calibration_config()

Return a named tuple consisting of
- `minibatch_size`: The size of each minibatch. Should be a multiple of the
  length of `sample_date_ranges`.
- `n_iterations`: The number of iterations.
- `rng`: RNG object for ensuring consistent calibration runs.
"""
function get_calibration_config()
    rng_seed = 42
    rng = Random.MersenneTwister(rng_seed)
    return (; minibatch_size = 2, n_iterations = 6, rng = rng)
end

"""
    get_ekp()

Return an `EnsembleKalmanProcess` object for calibration.
"""
function get_ekp()
    (; sample_date_ranges) = get_config()
    (; minibatch_size, rng) = get_calibration_config()
    # ensemble_size = 1
    obs_series = get_observation(sample_date_ranges, minibatch_size)
    # eki = EKP.EnsembleKalmanProcess(
    #     EKP.construct_initial_ensemble(rng, prior, ensemble_size),
    #     obs_series,
    #     EKP.TransformInversion(),
    #     verbose = true,
    #     scheduler = EKP.DataMisfitController(terminate_at = 100),
    # )

    # TODO: I feel like the constructor should be consistent, but they
    # are not :(
    eki = EKP.EnsembleKalmanProcess(obs_series, EKP.TransformUnscented(prior, impose_prior = true), verbose = true, rng = rng, scheduler = EKP.DataMisfitController(terminate_at = 100))
    ensemble_size = EKP.get_N_ens(eki)
    return eki
end

"""
    get_prior()


!!! warning "Check if parameters are overwritten!"
    Due to the difficulty of overwritting parameters in the land model, you must
    go to `model_interface.jl` and add which value is being overwritten and what
    prior to load it from. For example, if you add
    `EKP.constrained_gaussian("α_0", 0.64, 0.05, 0.2, 0.8)`, you should add
    ```julia
    α_0 = FT(get(calibrate_param_dict, "α_0", 0.64))
    ```
    and use `α_0` in the land model.
"""
function get_prior()
    priors = [
        # EKP.constrained_gaussian("Δα", 0.2, 0.1, 0.0, 1.0),
        # EKP.constrained_gaussian("k", 10, 5, 2, 25),
        # EKP.constrained_gaussian("beta_snow", 0.4, 0.2, 0.1, 0.8),
        # EKP.constrained_gaussian("x0_snow", 0.4, 0.2, 0.1, 0.8),
        EKP.constrained_gaussian("beta_snow_cover", 1.77, 0.2, 0.1, 2.5),
        EKP.constrained_gaussian("z0_snow_cover", 0.106, 0.03, 0.0, 0.3),
        EKP.constrained_gaussian("α_0", 0.7, 0.2, 0.0, 1.0),
        # EKP.constrained_gaussian("Δα", 0.7, 0.2, 0.0, 1.0);
        # EKP.constrained_gaussian("z0_snow", 0.106, 0.05, 0.01, 0.3),
    ]
    return EKP.combine_distributions(priors)
end

"""
    get_observation()

Get the obervation.

Before running this function, you should run `generate_observations.jl` to
make `land_observation_vector.jld2`.

!!! warning "Checking dates"
    There is no check for whether the `sample_date_ranges`
"""
function get_observation(sample_date_ranges, minibatch_size)
    observation_vector = JLD2.load_object(
        "experiments/better_calibration/land_observation_vector.jld2",
    )

    # TODO: Add a check for dates in case that the sample_date_ranges is modified, but
    # forgot to make a new observation
    # TODO: Might need to add a function in ClimaAnalysis and ClimaCalibrate to
    # make it easy to check for dates

    obs_series = EKP.ObservationSeries(
        Dict(
            "observations" => observation_vector,
            "names" => [
                string(Dates.year(start_date)) for
                (start_date, end_date) in sample_date_ranges
            ],
            "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
                length(observation_vector),
                minibatch_size,
            ),
        ),
    )
    return obs_series
end
