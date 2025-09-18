# Strong scaling runs for the soil paper

delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

include("strong_scaling_helper.jl")

update_drivers = false


run(snowy_land_integrator; info = false, update_drivers)
