tutorials = [
    "Running standalone component simulations" => [
        "Bucket" => "standalone/Bucket/bucket_tutorial.jl",
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
            "Changing soil parameterizations" => "standalone/Soil/changing_soil_parameterizations.jl",
        ],
        "Canopy" => [
            "Default Canopy" => "standalone/Canopy/default_canopy.jl",
            "Standalone Canopy" => "standalone/Canopy/canopy_tutorial.jl",
            "Changing canopy parameterizations" => "standalone/Canopy/changing_canopy_parameterizations.jl",
        ],
        "Snow" => [
            "Base tutorial" => "standalone/Snow/base_tutorial.jl",
            "Data tutorial" => "standalone/Snow/data_tutorial.jl",
        ],
    ],
    "Running Fluxnet simulations" => [
        "Canopy and soil" => "integrated/soil_canopy_fluxnet_tutorial.jl",
        "Canopy, soil, and snow" => "integrated/snowy_land_fluxnet_tutorial.jl",
        "Data processing" => "integrated/fluxnet_data.jl",
        "Visualization" => "integrated/fluxnet_vis.jl",
    ],
    "Running global simulations" => [
        "Bucket" => "global/bucket.jl",
        "Snow, soil, canopy" => "global/snowy_land.jl",
    ],
    "Calibrating a ClimaLand model" => [
        "Single site perfect model" => "calibration/minimal_working_example.jl",
        "Single site observations" => "calibration/minimal_working_example_obs.jl",
    ],
    "Running coupled simulations" =>
        ["Coupled bucket model" => "standalone/Bucket/coupled_bucket.jl"],
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
