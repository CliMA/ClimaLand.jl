#!/bin/bash

# implicit solver with flux_in=-1e-7, max_iters=1
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_16384 top_flux_-1e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_10000 top_flux_-1e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_12000 top_flux_-1e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_14000 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_8192 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_4096 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_2048 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_1024 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_512 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_256 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_128 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_64 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_16 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_4 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_1 top_flux_-1e-7

# implicit solver with flux_in=-1e-7, max_iters=2
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_16384 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_8192 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_4096 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_2048 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_1024 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_512 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_256 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_128 top_flux_-1e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_64 top_flux_-1e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_16 top_flux_-1e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_4 top_flux_-1e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_1 top_flux_-1e-7

# explicit solver with flux_in=-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_512 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_256 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_128 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_64 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_16 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_4 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_1 top_flux_-1e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_0p25 top_flux_-1e-7





# implicit solver with flux_in=-5e-7, max_iters=1
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_16384 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_8192 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_4096 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_2048 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_1024 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_512 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_256 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_128 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_64 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_16 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_4 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_1 top_flux_-5e-7

# implicit solver with flux_in=-5e-7, max_iters=2
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_16384 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_8192 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_4096 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_2048 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_1024 top_flux_-5e-7
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_512 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_256 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_128 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_64 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_16 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_4 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_2 dt_1 top_flux_-5e-7

# explicit solver with flux_in=-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_512 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_256 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_128 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_64 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_16 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_4 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_1 top_flux_-5e-7
# julia --project src/Soil/examples/imex_testing_cliargs.jl exp mod_picard iters_1 dt_0p25 top_flux_-5e-7
