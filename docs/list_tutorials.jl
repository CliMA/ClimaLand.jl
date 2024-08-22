tutorials = [
    "Running simulations" => [
        "Bucket LSM" => [
            "standalone/Bucket/bucket_tutorial.jl",
            "standalone/Bucket/coupled_bucket.jl",
        ],
        "Soil modeling" => [
            "Boundary conditions" => "standalone/Soil/boundary_conditions.jl",
            "Richards Equation" => "standalone/Soil/richards_equation.jl",
            "Energy and Hydrology" => "standalone/Soil/soil_energy_hydrology.jl",
            "Phase Changes" => "standalone/Soil/freezing_front.jl",
            "Layered Soil" => "standalone/Soil/layered_soil.jl",
            "Coarse Sand Evaporation" => "standalone/Soil/evaporation.jl",
            "Gilat Loess Evaporation" => "standalone/Soil/evaporation_gilat_loess.jl",
            "Bare soil site" => "standalone/Soil/sublimation.jl",
        ],
        "Canopy modeling" => [
            "Standalone Canopy" => "standalone/Canopy/canopy_tutorial.jl",
        ],
        "Integrated soil+canopy modeling" => [
            "Coupled Canopy and Soil" => "integrated/soil_canopy_tutorial.jl",
        ],
        "Snow Modeling" => [
            "standalone/Snow/base_tutorial.jl",
            "standalone/Snow/data_tutorial.jl",
        ],
    ],
    "For model developers" => [
        "Intro to standalone models" => "standalone/Usage/model_tutorial.jl",
        "Intro to multi-component models" => "standalone/Usage/LSM_single_column_tutorial.jl",
        "Intro to ClimaLand Domains" => "standalone/Usage/domain_tutorial.jl",
        "Intro to forced site-level runs" => "shared_utilities/driver_tutorial.jl",
    ],
]
