tutorials = [
    "Running Fluxnet simulations" => [
        "Canopy and soil" => "integrated/soil_canopy_fluxnet_tutorial.jl",
        "Canopy, soil, and snow" => "integrated/snowy_land_fluxnet_tutorial.jl",
        "Data processing" =>
            "Fluxnet forcing and comparison data" => "integrated/fluxnet_data.jl",
        "Visualization" =>
            "Fluxnet simulation visualization" => "integrated/fluxnet_vis.jl",
    ],
    "Running global simulations" => [
        "Bucket land model" => [
            "standalone/Bucket/bucket_tutorial.jl",
            "standalone/Bucket/coupled_bucket.jl",
        ],
    ],
    "Running standalone component simulations" => [
        "Soil" => [
            "Boundary conditions" => "standalone/Soil/boundary_conditions.jl",
            "Richards Equation" => "standalone/Soil/richards_equation.jl",
            "Energy and Hydrology" => "standalone/Soil/soil_energy_hydrology.jl",
            "Phase Changes" => [
                "standalone/Soil/freezing_front.jl",
                "standalone/Soil/phase_change_analytic.jl",
            ],
            "Layered Soil" => "standalone/Soil/layered_soil.jl",
            "Coarse Sand Evaporation" => "standalone/Soil/evaporation.jl",
            "Gilat Loess Evaporation" => "standalone/Soil/evaporation_gilat_loess.jl",
            "Bare soil site" => "standalone/Soil/sublimation.jl",
        ],
        "Canopy" => [
            "Standalone Canopy" => "standalone/Canopy/canopy_tutorial.jl",
        ],
        "Snow" => [
            "standalone/Snow/base_tutorial.jl",
            "standalone/Snow/data_tutorial.jl",
        ],
    ],
    "Calibrating a ClimaLand model" => [
        "Single site perfect model" => "calibration/minimal_working_example.jl",
        "Single site observations" => "calibration/minimal_working_example_obs.jl",
    ],
    "For model developers" => [
        "Intro to standalone models" => "standalone/Usage/model_tutorial.jl",
        "Intro to multi-component models" => [
            "Single column tutorial" => "standalone/Usage/LSM_single_column_tutorial.jl",
            "Adjusting boundary conditions for the soil" => "integrated/handling_soil_fluxes.jl",
            "Adjusting boundary conditions for the snow" => "integrated/handling_snow_fluxes.jl",
        ],
        "Intro to ClimaLand Domains" => "standalone/Usage/domain_tutorial.jl",
        "Intro to forced site-level runs" => "shared_utilities/driver_tutorial.jl",
        "Intro to implicit/explicit timestepping" => "shared_utilities/timestepping.jl",
    ],
]
