# Runs for the soil paper

delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

include("weak_scaling_helper.jl")

update_drivers = false



run(
    snow_integrator;
    label = "Snow",
    plot_attrs = "color=black, mark=+",
    info = false,
    update_drivers,
    tf = 86400.0 * 2,
)

run(
    bucket_integrator;
    label = "Bucket",
    plot_attrs = "color=red, mark=square*, mark size = 1.8",
    info = false,
    update_drivers,
    tf = 86400.0 * 2,
)

run(
    richards_integrator;
    label = "Richards",
    plot_attrs = "color=green!50!black, mark=triangle*, mark size = 2.5",
    info = false,
    update_drivers,
    tf = 86400.0 * 2,
)



run(
    canopy_integrator;
    label = "Canopy",
    plot_attrs = "color=purple!80!black, mark=diamond*, mark size=2.5",
    info = false,
    update_drivers,
)

run(
    soil_integrator;
    label = "Soil",
    plot_attrs = "color=red!40!white, mark=halfcircle*",
    info = false,
    update_drivers,
)

run(
    soil_canopy_integrator;
    label = "Soil-Canopy",
    plot_attrs = "color=orange!80!black, mark=pentagon*",
    info = false,
    update_drivers,
)

run(
    snowy_land_integrator;
    label = "Soil-Canopy-Snow",
    plot_attrs = "color=blue, mark=*",
    info = false,
    update_drivers,
)
